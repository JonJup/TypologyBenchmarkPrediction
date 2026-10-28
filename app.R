### PULSE benchmark predictor: quantile random forest app ###
#
# Serves the QRF benchmark models written by 11_fit_qrf.R, unchanged. The model
# files sit in a folder next to this script (MODEL_DIR, default "models"):
#   <taxon>__<metric>.rds            list(model = ranger QRF, recipe, best, predictors)
#   <taxon>__<metric>_intervals.rds  optional; its held-out PIT values enable calibration
# Taxa and metrics are discovered from the file names at startup.
# The prepped recipe retains the training set (recipe$template: predictors, value,
# data.set). The app uses it as the reference typologies for imputation, and to
# trace each tree's prediction (one training value per tree) to its data set.
#
# Uncertainty reported by the app
#   Data sets + trees : cluster bootstrap of the per-tree predictions; trees and the
#                       training data sets that supply their values are resampled.
#                       The forest structure is held fixed.
#   Calibration       : (with _intervals.rds) quantiles and percentiles mapped
#                       through the empirical CDF of held-out PIT values; held-out
#                       data sets are resampled in the same bootstrap.
#   Imputation        : unspecified inputs are varied over their plausible range
#                       among the training typologies, one at a time and jointly.
#
# Reduced forests (12_reduce_qrf.R) have the same format plus $reduction; the app
# reads them unchanged and shows the number of trees kept in Model Info.
#
# Memory. A loaded model occupies its forest; each prediction adds a temporary
# working copy of about twice the forest. At most QRF_CACHE models are held in
# memory, across all sessions: a session keeps only the model's metadata and takes
# the forest from the shared cache when it predicts.
#
# Environment variables
#   MODEL_DIR        folder with the model files (default "models")
#   QRF_CACHE        models kept in memory (default 1)
#   QRF_THREADS      ranger threads (default 2)
#   FULL_MODELS_URL  link to the full models, shown in Model Info (optional)


# Startup diagnostics -----------------------------------------------------

cat("\n\nQRF Prediction App - Starting...\n\n")

# recipes is not used directly, but the model files contain a prepped recipe,
# which readRDS() can only restore with the package installed.
required_packages <- c("shiny", "bslib", "ranger", "recipes", "data.table", "DT", "ggplot2")

missing_packages <- character(0)
for (pkg in required_packages) {
        ok <- requireNamespace(pkg, quietly = TRUE)
        cat(if (ok) "OK:     " else "MISSING:", pkg, "\n")
        if (!ok) missing_packages <- c(missing_packages, pkg)
}
if (length(missing_packages) > 0) {
        stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
             "\nInstall with: install.packages(c('", paste(missing_packages, collapse = "', '"), "'))")
}
cat("\nAll packages available. Loading...\n\n")

suppressPackageStartupMessages({
        library(shiny)
        library(bslib)
        library(ranger)
        library(data.table)
        library(DT)
        library(ggplot2)
})

if (!exists("%||%", baseenv())) `%||%` <- function(x, y) if (is.null(x)) y else x



# Configuration -----------------------------------------------------------

# Folder with the model files, relative to this script (or an absolute path).
MODEL_DIR <- Sys.getenv("MODEL_DIR", "models")
if (!grepl("^(/|~|[A-Za-z]:)", MODEL_DIR)) MODEL_DIR <- file.path(getwd(), MODEL_DIR)
MODEL_DIR <- normalizePath(MODEL_DIR, mustWork = FALSE)

env_int <- function(name, default) {
        v <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
        if (is.na(v) || v < 1L) default else v
}
N_THREADS       <- env_int("QRF_THREADS", 2L)
MODEL_CACHE_MAX <- env_int("QRF_CACHE", 1L)   # models kept in memory, across all sessions
FULL_MODELS_URL <- Sys.getenv("FULL_MODELS_URL", "")
REF_MAX_ROWS    <- 50000L            # training typologies kept for KNN imputation

# ---- Uncertainty ----------------------------------------------------------------
CI_LEVEL         <- 0.90             # level of all reported uncertainty intervals
RANGE_PROBS      <- c(0.10, 0.90)    # plausible range of an imputed input among reference typologies
N_REF_DRAWS      <- 200L             # reference typologies drawn to propagate imputation jointly
N_BOOT           <- 1000L            # bootstrap replicates (trees, training and held-out data sets)
PIT_BINS         <- 1000L            # resolution of the calibration map
CALIBRATION_SETS <- c("Test set", "Grouped CV (out-of-fold)")   # held-out PIT used for calibration
MIN_CAL_GROUPS   <- 5L               # fewer held-out data sets -> no calibration
LOW_SUPPORT      <- 3                # effective number of data sets below which support is flagged
MAX_UNMATCHED    <- 0.05             # tolerated share of tree values not found among training values
# Imputation impact: joint spread of the median (RANGE_PROBS) relative to the
# width of the benchmark's own RANGE_PROBS interval.
IMPACT_BREAKS    <- c(negligible = 0.10, minor = 0.33, moderate = 0.66)
N_TORNADO        <- 10L              # imputed inputs shown in the effect plot


# Model discovery ---------------------------------------------------------
# Every <taxon>__<metric>.rds in MODEL_DIR is a model, except the other outputs
# of 11_fit_qrf.R (save_out() suffixes), which may sit in the same folder.
AUX_SUFFIXES <- c("ale", "datasets", "intervals", "lofo", "metrics",
                  "quantile_importance", "run_info", "split", "tuning")

# Evaluation metrics: display names, in menu order. Keys are the metric part of
# the file names (clean_name() in 11_fit_qrf.R turns spaces into "_"). Metrics
# found in the folder but not listed here are shown after these, by file name.
METRIC_LABELS <- c(
        Classification_Strength = "Classification Strength",
        ANOSIM_R                = "ANOSIM R",
        AucZeta                 = "AucZeta",
        PERMANOVA_R2            = "PERMANOVA R2",
        PERMANOVA_Fuzzy_R2      = "PERMANOVA Fuzzy R2",
        Fuzzy_Mantel            = "Fuzzy Mantel"
)
metric_display_label <- function(m) {
        vapply(m, function(x) if (x %in% names(METRIC_LABELS)) METRIC_LABELS[[x]] else gsub("_", " ", x), "",
               USE.NAMES = FALSE)
}

discover_models <- function(dir) {
        stems <- sub("\\.rds$", "", list.files(dir, pattern = "__.+\\.rds$"))
        stems <- stems[!grepl(paste0("_(", paste(AUX_SUFFIXES, collapse = "|"), ")$"), stems)]
        if (!length(stems)) return(data.table(taxon = character(), metric = character()))
        # clean_name() in 11_fit_qrf.R never yields "__" inside a taxon or metric
        data.table(taxon = sub("__.*$", "", stems), metric = sub("^.*?__", "", stems))[order(taxon, metric)]
}

MODELS <- discover_models(MODEL_DIR)
cat("Model folder:", MODEL_DIR, "\n")
cat("Models kept in memory:", MODEL_CACHE_MAX, "| ranger threads:", N_THREADS, "\n")
if (nrow(MODELS)) {
        for (t in unique(MODELS$taxon)) cat(sprintf("  %-20s %s\n", t, paste(MODELS[taxon == t, metric], collapse = ", ")))
        unknown <- setdiff(MODELS$metric, names(METRIC_LABELS))
        if (length(unknown)) cat("  Note: metrics without a display name in METRIC_LABELS:", paste(unknown, collapse = ", "), "\n")
} else {
        cat("  WARNING: no model files (<taxon>__<metric>.rds) found.\n")
}
cat("\n")

# Display names and icons for known taxa; others are title-cased.
KNOWN_TAXA <- list(
        diatoms            = list(label = "Diatoms",            icon = "leaf"),
        fishes             = list(label = "Fishes",             icon = "fish"),
        fish               = list(label = "Fishes",             icon = "fish"),
        macroinvertebrates = list(label = "Macroinvertebrates", icon = "bug"),
        macrophytes        = list(label = "Macrophytes",        icon = "seedling")
)
taxon_info <- function(t) {
        KNOWN_TAXA[[tolower(t)]] %||% list(label = tools::toTitleCase(gsub("_", " ", t)), icon = "water")
}

# Metrics available for a taxon, as c(label = file metric) for selectInput()
get_available_metrics <- function(taxon_group = NULL) {
        if (is.null(taxon_group) || !nzchar(taxon_group)) return(character(0))
        m <- MODELS[taxon == taxon_group, metric]
        m <- m[order(match(m, names(METRIC_LABELS)), m)]
        setNames(m, metric_display_label(m))
}


# Predictor specification -------------------------------------------------
# The predictors of 11_fit_qrf.R (predictor_vec). Defaults are only used when a
# model file carries no training data; otherwise unspecified inputs are imputed
# from the training typologies.
PREDICTOR_SPEC <- list(
        # Typology
        n_types = list(label = "Number of types", group = "Typology", is_integer = TRUE,
                       min = 2, max = 50, value = 5, step = 1,
                       help = "Number of types in the typology."),
        n_vars = list(label = "Number of typology variables", group = "Typology", is_integer = TRUE,
                      min = 1, max = 50, value = 5, step = 1,
                      help = "Number of environmental variables the typology is built on."),
        type_evenness = list(label = "Type evenness", group = "Typology",
                             min = 0, max = 1, value = 0.9, step = 0.01,
                             help = "Evenness of the type sizes (sites per type); 1 = all types equally large."),
        env_r2_after = list(label = "Environmental R\u00b2", group = "Typology",
                            min = 0, max = 1, value = 0.3, step = 0.01,
                            help = "Share of the variance of all standardised environmental variables explained by the types (all variables, not only the typology variables)."),
        env_ch_after = list(label = "Calinski\u2013Harabasz index", group = "Typology",
                            min = 0, max = 1000, value = 20, step = 1,
                            help = "Calinski\u2013Harabasz index of the types, in the space of all standardised environmental variables."),
        asw_after = list(label = "Average silhouette width", group = "Typology",
                         min = -1, max = 1, value = 0.15, step = 0.01,
                         help = "Average silhouette width of the types, in the space of all standardised environmental variables."),
        
        # Spatial
        samples = list(label = "Number of sites", group = "Spatial", is_integer = TRUE,
                       min = 3, max = 10000, value = 100, step = 1,
                       help = "Number of sites in the typology."),
        min_distance = list(label = "Minimum distance (km)", group = "Spatial",
                            min = 0, max = 10000, value = 1, step = 0.1,
                            help = "Smallest pairwise distance between sites."),
        median_distance = list(label = "Median distance (km)", group = "Spatial",
                               min = 0, max = 10000, value = 100, step = 1,
                               help = "Median pairwise distance between sites."),
        
        # Community (variation partitioning of the observed community)
        n_taxa_used = list(label = "Number of taxa", group = "Community", is_integer = TRUE,
                           min = 2, max = 5000, value = 50, step = 1,
                           help = "Number of taxa that occur at least once at the sites."),
        vp_env_pure = list(label = "Pure environmental fraction [a]", group = "Community",
                           min = -1, max = 1, value = 0.05, step = 0.01,
                           help = "Adjusted R\u00b2 of the typology variables after accounting for space (variation partitioning of Jaccard dissimilarities with the selected MEMs)."),
        vp_shared = list(label = "Shared fraction [b]", group = "Community",
                         min = -1, max = 1, value = 0.03, step = 0.01,
                         help = "Adjusted R\u00b2 shared by the typology variables and the selected MEMs."),
        vp_space_pure = list(label = "Pure spatial fraction [c]", group = "Community",
                             min = -1, max = 1, value = 0.03, step = 0.01,
                             help = "Adjusted R\u00b2 of the selected MEMs after accounting for the typology variables."),
        
        # Taxonomic resolution
        species_rank = list(label = "Species rank (share)", group = "Taxonomic",
                            min = 0, max = 1, value = 0.7, step = 0.01,
                            help = "Share of taxa identified to species."),
        genus_rank = list(label = "Genus rank (share)", group = "Taxonomic",
                          min = 0, max = 1, value = 0.2, step = 0.01,
                          help = "Share of taxa identified to genus."),
        family_rank = list(label = "Family rank (share)", group = "Taxonomic",
                           min = 0, max = 1, value = 0.08, step = 0.01,
                           help = "Share of taxa identified to family."),
        higher_rank = list(label = "Higher rank (share)", group = "Taxonomic",
                           min = 0, max = 1, value = 0.02, step = 0.01,
                           help = "Share of taxa identified above family."),
        
        # Sampling dates
        mean_doy = list(label = "Mean sampling day of year", group = "Sampling",
                        min = 1, max = 366, value = 180, step = 1,
                        help = "Mean sampling date, as day of the year."),
        circ_sd_days = list(label = "Spread of sampling dates (days)", group = "Sampling",
                            min = 0, max = 183, value = 30, step = 1,
                            help = "Circular standard deviation of the sampling dates, in days.")
)
PREDICTOR_GROUPS <- unique(vapply(PREDICTOR_SPEC, `[[`, "", "group"))

spec_label   <- function(v) vapply(v, function(x) PREDICTOR_SPEC[[x]]$label %||% x, "")
spec_default <- function(v) vapply(v, function(x) PREDICTOR_SPEC[[x]]$value %||% 0, 0)
is_int_var   <- function(v) vapply(v, function(x) isTRUE(PREDICTOR_SPEC[[x]]$is_integer), TRUE)

# Constraints on imputed values. Specified values are never changed.
#   equal : the group sums to `total`
#   max   : the group sums to at most `total` (adjusted R2 fractions; [d] = 1 - a - b - c)
SUM_CONSTRAINTS <- list(
        taxonomic = list(vars = c("species_rank", "genus_rank", "family_rank", "higher_rank"),
                         type = "equal", total = 1),
        variation = list(vars = c("vp_env_pure", "vp_shared", "vp_space_pure"),
                         type = "max", total = 1)
)
# Ordered chains: each element <= the next
ORDER_CONSTRAINTS <- list(distance = c("min_distance", "median_distance"))


# Formatting helpers ------------------------------------------------------

fmt <- function(x, digits = 3) {
        ifelse(is.na(x), "NA",
               ifelse(abs(x) >= 100, formatC(x, format = "f", digits = 0),
                      formatC(x, format = "f", digits = digits)))
}
fmt_in <- function(x) trimws(formatC(x, digits = 3, format = "fg"))
ordinal <- function(n) {
        n <- as.integer(n)
        suffix <- if (n %% 100 %in% 11:13) "th" else switch(as.character(n %% 10),
                                                            "1" = "st", "2" = "nd", "3" = "rd", "th")
        paste0(n, suffix)
}
pct_label <- function(p) {
        x <- 100 * p
        if (x < 1) "<1st" else if (x > 99) ">99th" else ordinal(round(x))
}
pct_range <- function(ci) sprintf("%d\u2013%d", as.integer(round(100 * ci[1])), as.integer(round(100 * ci[2])))
level_pct <- function(level) as.integer(round(100 * level))


# Model memo --------------------------------------------------------------

.memo <- new.env(parent = emptyenv())
.memo$items <- list()
memo_set <- function(key, value) {
        items <- c(setNames(list(value), key), .memo$items[setdiff(names(.memo$items), key)])
        .memo$items <- head(items, MODEL_CACHE_MAX)
        invisible(value)
}
# Drop the least recently used models until n remain.
memo_trim <- function(n) {
        if (length(.memo$items) > n) {
                .memo$items <- head(.memo$items, max(0L, n))
                invisible(gc(verbose = FALSE))
        }
}


# Load model, training data and (optional) held-out PIT; memoised across sessions.
load_benchmark <- function(taxon, metric) {
        key <- paste0(taxon, "__", metric)
        if (!is.null(.memo$items[[key]])) return(memo_set(key, .memo$items[[key]]))
        
        # Make room first, so that loading never holds more than MODEL_CACHE_MAX models
        memo_trim(MODEL_CACHE_MAX - 1L)
        path <- file.path(MODEL_DIR, paste0(key, ".rds"))
        if (!file.exists(path)) stop("Model file not found: ", path)
        cat("  Reading:", path, sprintf("(%.0f MB)", file.size(path) / 1e6), "\n")
        obj <- readRDS(path)
        if (!is.list(obj) || !inherits(obj$model, "ranger")) {
                stop(key, " is not a model file written by 11_fit_qrf.R (a list with a ranger $model).")
        }
        if (is.null(obj$model$random.node.values)) stop(key, " was not fitted with quantreg = TRUE.")
        preds <- obj$model$forest$independent.variable.names
        
        training <- tryCatch(prepare_training(obj$recipe, preds), error = function(e) {
                cat("  Training data unusable:", conditionMessage(e), "\n"); NULL })
        model     <- obj$model
        reduction <- obj$reduction                  # NULL for full forests
        rm(obj)
        
        int_path <- file.path(MODEL_DIR, paste0(key, "_intervals.rds"))
        ints <- if (file.exists(int_path)) tryCatch(readRDS(int_path), error = function(e) NULL) else NULL
        calibration <- tryCatch(build_calibration(ints$pit), error = function(e) NULL)
        coverage <- if (!is.null(ints$test_intervals)) as.data.table(ints$test_intervals) else NULL
        rm(ints); invisible(gc(verbose = FALSE))
        
        memo_set(key, list(key = key, taxon = taxon, metric = metric, model = model, predictors = preds,
                           n_trees = model$num.trees, training = training,
                           calibration = calibration, coverage = coverage, reduction = reduction))
}

# Session copy without the forest and the reference typologies (the large parts).
benchmark_light <- function(bm) {
        if (is.null(bm)) return(NULL)
        bm$model <- NULL
        if (!is.null(bm$training)) {
                bm$training$reference <- NULL
                bm$training$donors    <- NULL
        }
        bm
}

# Training set from the prepped recipe (prep(retain = TRUE) keeps it as $template).
prepare_training <- function(rec, preds) {
        tmpl <- rec$template
        if (is.null(tmpl)) return(NULL)
        need <- c(preds, "value", "data.set")
        if (!all(need %in% names(tmpl))) return(NULL)
        tmpl <- as.data.frame(tmpl)[, need]
        tmpl <- tmpl[stats::complete.cases(tmpl), , drop = FALSE]
        n    <- nrow(tmpl)
        keep <- if (n > REF_MAX_ROWS) with_seed(1L, sort(sample.int(n, REF_MAX_ROWS))) else seq_len(n)
        X    <- tmpl[, preds, drop = FALSE]
        list(reference  = X[keep, , drop = FALSE],
             medians    = vapply(X, median, 0),
             min        = vapply(X, min, 0),
             max        = vapply(X, max, 0),
             n_rows     = n,
             n_datasets = length(unique(tmpl$data.set)),
             donors     = build_donor_lookup(tmpl$value, tmpl$data.set))
}


# Data-set support -----------------------------------------------------------
# Each tree predicts one training value (ranger's random.node.values). Matching it
# to the training values identifies the data set it came from. Where a value
# occurs in several data sets (ties, e.g. count metrics), one is drawn in
# proportion to its frequency.

build_donor_lookup <- function(values, groups) {
        f  <- factor(groups)
        dt <- data.table(value = values, g = as.integer(f))[, .N, by = .(value, g)]
        setorder(dt, value, -N)
        dt[, `:=`(first = .I[1L], n_g = .N), by = value]
        u <- unique(dt[, .(value, first, n_g)])
        list(values = u$value, first = u$first, n_g = u$n_g, g = dt$g, N = dt$N,
             groups = levels(f), n_groups = nlevels(f))
}

donor_groups <- function(lk, v, seed = 1L) {
        i <- match(v, lk$values)
        g <- lk$g[lk$first[i]]
        amb <- which(!is.na(i) & lk$n_g[i] > 1L)
        if (length(amb)) with_seed(seed, for (t in amb) {
                pos <- lk$first[i[t]] + seq_len(lk$n_g[i[t]]) - 1L
                g[t] <- lk$g[pos][sample.int(length(pos), 1L, prob = lk$N[pos])]
        })
        list(g = g, ambiguous = length(amb) / length(v), unmatched = mean(is.na(i)))
}

support_stats <- function(dg, n_groups) {
        cnt <- tabulate(dg$g[!is.na(dg$g)], n_groups)
        s   <- cnt[cnt > 0] / sum(cnt)
        list(n_datasets = length(s), n_eff = 1 / sum(s^2), top_share = max(s),
             ambiguous = dg$ambiguous, unmatched = dg$unmatched)
}


# Calibration -------------------------------------------------------------
# The held-out PIT values are binned per data set (H: data sets x bins). The
# calibration map C(p) is the share of held-out typologies whose PIT is <= p,
# i.e. the share falling below the forest's p-quantile. Bootstrap maps resample
# held-out data sets (rows of H) with replacement.

cum_share <- function(counts) c(0, cumsum(counts)) / sum(counts)

build_calibration <- function(pit) {
        if (is.null(pit)) return(NULL)
        pit <- as.data.table(pit)
        if (!all(c("set", "data.set", "pit") %in% names(pit))) return(NULL)
        pit <- pit[set %in% CALIBRATION_SETS & is.finite(pit)]
        if (uniqueN(pit$data.set) < MIN_CAL_GROUPS) return(NULL)
        
        edges <- seq(0, 1, length.out = PIT_BINS + 1L)
        bin   <- findInterval(pit$pit, edges, rightmost.closed = TRUE, all.inside = TRUE)
        H     <- unclass(table(factor(pit$data.set), factor(bin, levels = seq_len(PIT_BINS))))
        storage.mode(H) <- "double"
        
        G <- nrow(H)
        W <- with_seed(1L, vapply(seq_len(N_BOOT), function(b)
                as.numeric(tabulate(sample.int(G, G, replace = TRUE), G)), numeric(G)))   # G x B
        boot <- t(apply(crossprod(W, H), 1, cum_share))                                   # B x (bins + 1)
        
        list(edges = edges, cum = cum_share(colSums(H)), boot = boot,
             n_datasets = G, n_rows = nrow(pit), sets = unique(pit$set))
}

# Share of held-out typologies below the forest's p-quantile
cal_cdf <- function(cum, edges, p) approx(edges, cum, xout = p, rule = 2, ties = "ordered")$y

# Forest level whose quantile has a share tau of held-out typologies below it
cal_inv <- function(cum, edges, tau) {
        j  <- pmin(length(cum), pmax(2L, findInterval(tau, cum, left.open = TRUE) + 1L))
        lo <- cum[j - 1L]; hi <- cum[j]
        w  <- ifelse(hi > lo, (tau - lo) / (hi - lo), 0.5)
        edges[j - 1L] + pmin(1, pmax(0, w)) * (edges[j] - edges[j - 1L])
}


# Forest helpers ----------------------------------------------------------

# Evaluate expr with a fixed seed and restore the caller's RNG state.
with_seed <- function(seed, expr) {
        old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv()) else NULL
        on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
        set.seed(seed)
        expr
}

# Per-tree predictions (n x trees): one training value per tree, drawn from the
# terminal node at fit time. Their empirical distribution is the QRF predictive
# distribution; quantile(type = 7) of a row reproduces predict(type = "quantiles").
tree_values <- function(model, X) {
        p <- predict(model, data = X, type = "quantiles", what = function(v) v,
                     num.threads = N_THREADS)$predictions
        p <- as.matrix(p)
        if (nrow(p) != nrow(X)) p <- t(p)
        p
}

# Type-7 quantiles and mid-rank percentile of a sorted vector
q7s <- function(vs, p) {
        n <- length(vs); p <- pmin(1, pmax(0, p))
        h <- (n - 1) * p + 1; lo <- floor(h); hi <- pmin(n, lo + 1)
        vs[lo] + (h - lo) * (vs[hi] - vs[lo])
}
mid_pct <- function(vs, x) {
        below <- findInterval(x, vs, left.open = TRUE)    # values < x
        upto  <- findInterval(x, vs)                      # values <= x
        (below + 0.5 * (upto - below)) / length(vs)
}

# Quantiles and observed percentile with uncertainty, in one bootstrap loop.
# Replicate b resamples the trees and (if groups are known) the training data
# sets; a tree's value counts (tree multiplicity) x (multiplicity of its data
# set). With calibration, replicate b also uses the b-th bootstrap calibration map.
summarise_benchmark <- function(v, taus, obs = NULL, cal = NULL, groups = NULL, n_groups = NULL, seed = 1L) {
        vs <- sort(v); n <- length(vs); a <- (1 - CI_LEVEL) / 2; L <- length(taus)
        has_obs <- !is.null(obs); has_cal <- !is.null(cal); has_grp <- !is.null(groups)
        ord <- order(v); v_o <- v[ord]
        g_o <- if (has_grp) groups[ord] else NULL
        
        boot <- with_seed(seed, lapply(seq_len(N_BOOT), function(b) {
                w <- tabulate(sample.int(n, n, replace = TRUE), n)
                if (has_grp) {
                        m <- tabulate(sample.int(n_groups, n_groups, replace = TRUE), n_groups)
                        w <- w * ifelse(is.na(g_o), 1L, m[g_o])
                }
                if (!any(w > 0)) return(list(q = rep(NA_real_, L), qc = rep(NA_real_, L), p = NA_real_, pc = NA_real_))
                vb  <- rep.int(v_o, w)                     # sorted, because v_o is sorted
                out <- list(q = q7s(vb, taus))
                if (has_obs) out$p <- mid_pct(vb, obs)
                if (has_cal) {
                        cb <- cal$boot[b, ]
                        out$qc <- q7s(vb, cal_inv(cb, cal$edges, taus))
                        if (has_obs) out$pc <- cal_cdf(cb, cal$edges, out$p)
                }
                out
        }))
        ci <- function(nm) {
                M <- matrix(unlist(lapply(boot, `[[`, nm)), ncol = N_BOOT)
                t(apply(M, 1, quantile, probs = c(a, 1 - a), names = FALSE, na.rm = TRUE))
        }
        
        qt <- data.table(tau = taus, value = q7s(vs, taus))
        qt[, c("lo", "hi") := as.data.table(ci("q"))]
        if (has_cal) {
                qt[, cal_value := q7s(vs, cal_inv(cal$cum, cal$edges, taus))]
                qt[, c("cal_lo", "cal_hi") := as.data.table(ci("qc"))]
        }
        
        observed <- NULL
        if (has_obs) {
                p <- mid_pct(vs, obs)
                observed <- list(value = obs, pct = p, pct_ci = ci("p")[1, ])
                if (has_cal) {
                        observed$cal_pct <- cal_cdf(cal$cum, cal$edges, p)
                        observed$cal_ci  <- ci("pc")[1, ]
                }
        }
        list(quantiles = qt, observed = observed, median = q7s(vs, 0.5),
             range = q7s(vs, RANGE_PROBS), n_trees = n, calibrated = has_cal, grouped = has_grp)
}

# Tree values -> data-set groups for the bootstrap, plus support statistics
tree_groups <- function(bm, v) {
        lk <- bm$training$donors
        if (is.null(lk)) return(NULL)
        dg <- donor_groups(lk, v)
        st <- support_stats(dg, lk$n_groups)
        list(groups = if (dg$unmatched <= MAX_UNMATCHED) dg$g else NULL, n_groups = lk$n_groups, support = st)
}

# Specified inputs outside the training range (QRF predictions are flat there)
outside_range <- function(bm, specified_vals) {
        tr <- bm$training
        if (is.null(tr) || !length(specified_vals)) return(NULL)
        v  <- intersect(names(specified_vals), names(tr$min))
        x  <- as.numeric(unlist(specified_vals[v]))
        bad <- x < tr$min[v] | x > tr$max[v]
        if (!any(bad)) return(NULL)
        data.table(var = v[bad], value = x[bad], min = tr$min[v][bad], max = tr$max[v][bad])
}

density_from_values <- function(v, from = min(v), to = max(v)) {
        if (length(v) < 3 || diff(range(v)) <= 0) return(NULL)
        bw <- 1.2 * bw.nrd0(v)
        if (!is.finite(bw) || bw <= 0) bw <- diff(range(v)) / 50
        d <- density(v, bw = bw, n = 512, from = from, to = to)
        data.frame(x = d$x, density = d$y)
}


# Imputation --------------------------------------------------------------

# Enforce SUM_CONSTRAINTS, ORDER_CONSTRAINTS and integer rounding on all
# columns not in `fixed` (row-wise; df holds one row per input combination).
enforce_constraints <- function(df, fixed) {
        for (g in SUM_CONSTRAINTS) {
                vars <- intersect(g$vars, names(df))
                free <- setdiff(vars, fixed)
                if (!length(free)) next
                fix       <- intersect(vars, fixed)
                fixed_sum <- if (length(fix)) rowSums(df[, fix, drop = FALSE]) else rep(0, nrow(df))
                remaining <- pmax(0, g$total - fixed_sum)
                M         <- as.matrix(df[, free, drop = FALSE])
                s_free    <- rowSums(M)
                if (g$type == "equal") {
                        pos <- s_free > 0
                        M[pos, ]  <- M[pos, , drop = FALSE] * (remaining[pos] / s_free[pos])
                        M[!pos, ] <- remaining[!pos] / length(free)
                } else {
                        over <- s_free > remaining & s_free > 0
                        M[over, ] <- M[over, , drop = FALSE] * (remaining[over] / s_free[over])
                }
                df[, free] <- M
        }
        for (chain in ORDER_CONSTRAINTS) {
                chain <- intersect(chain, names(df))
                if (length(chain) < 2) next
                for (i in seq(2, length(chain))) if (!chain[i] %in% fixed)
                        df[[chain[i]]] <- pmax(df[[chain[i]]], df[[chain[i - 1]]])
                for (i in rev(seq_len(length(chain) - 1))) if (!chain[i] %in% fixed)
                        df[[chain[i]]] <- pmin(df[[chain[i]]], df[[chain[i + 1]]])
        }
        for (v in setdiff(names(df), fixed)) if (is_int_var(v)) df[[v]] <- round(df[[v]])
        df
}

# Keep a value of `v` inside the bounds set by specified members of its chain.
clamp_order <- function(v, x, row, specified) {
        for (chain in ORDER_CONSTRAINTS) {
                i <- match(v, chain)
                if (is.na(i)) next
                lower <- intersect(chain[seq_len(i - 1)], specified)
                upper <- intersect(chain[-seq_len(i)], specified)
                if (length(lower)) x <- pmax(x, max(unlist(row[1, lower])))
                if (length(upper)) x <- pmin(x, min(unlist(row[1, upper])))
        }
        x
}

# k nearest reference typologies in the standardised specified inputs
knn_prep <- function(ref, vars) {
        M  <- as.matrix(ref[, vars, drop = FALSE])
        mu <- colMeans(M)
        s  <- apply(M, 2, sd); s[!is.finite(s) | s == 0] <- 1
        list(vars = vars, mu = mu, s = s, Z = sweep(sweep(M, 2, mu, "-"), 2, s, "/"))
}
knn_index <- function(prep, x, k) {
        q  <- (unname(x[prep$vars]) - prep$mu) / prep$s
        d2 <- rowSums(sweep(prep$Z, 2, q, "-")^2)
        order(d2)[seq_len(min(k, length(d2)))]
}

# Resolve one input combination. specified_vals: named list or vector of the
# inputs the user gave. Returns the full input row, the imputed inputs and the
# pool of reference typologies the imputed values come from.
resolve_inputs <- function(bm, specified_vals, method = "conditional", k = 50, knn_cache = NULL) {
        preds     <- bm$predictors
        specified <- intersect(names(specified_vals), preds)
        imputed   <- setdiff(preds, specified)
        ref       <- bm$training$reference
        vals      <- setNames(rep(NA_real_, length(preds)), preds)
        vals[specified] <- as.numeric(unlist(specified_vals[specified]))
        pool <- NULL
        method_used <- method
        
        if (length(imputed)) {
                if (is.null(ref)) {
                        method_used <- "defaults"
                        vals[imputed] <- spec_default(imputed)
                } else if (method == "conditional" && length(specified)) {
                        key  <- paste(sort(specified), collapse = "|")
                        prep <- if (!is.null(knn_cache)) knn_cache[[key]] else NULL
                        if (is.null(prep)) {
                                prep <- knn_prep(ref, specified)
                                if (!is.null(knn_cache)) knn_cache[[key]] <- prep
                        }
                        pool <- ref[knn_index(prep, vals, k), , drop = FALSE]
                        vals[imputed] <- colMeans(pool[, imputed, drop = FALSE])
                } else {
                        method_used <- "marginal"        # also: KNN with nothing to condition on
                        pool <- ref
                        vals[imputed] <- bm$training$medians[imputed]
                }
        }
        row <- enforce_constraints(as.data.frame(as.list(vals), check.names = FALSE), fixed = specified)
        list(row = row, specified = specified, imputed = imputed, pool = pool, method = method_used)
}

# Effect of the imputed inputs on the benchmark.
#   One at a time: each imputed input at the RANGE_PROBS quantiles of the pool,
#                  the others at their imputed values.
#   Jointly:       up to N_REF_DRAWS reference typologies from the pool supply all
#                  imputed inputs together (their covariance is kept). The pooled
#                  per-tree predictions form the benchmark that integrates over
#                  the imputed inputs.
imputation_analysis <- function(bm, inp, v_base, taus, obs = NULL, seed = 1L) {
        pool <- inp$pool; imputed <- inp$imputed; specified <- inp$specified; base <- inp$row
        if (is.null(pool) || !length(imputed)) return(NULL)
        
        rng <- t(vapply(imputed, function(x) {
                r <- q7s(sort(pool[[x]]), RANGE_PROBS)
                if (is_int_var(x)) r <- round(r)
                clamp_order(x, r, base, specified)
        }, numeric(2)))
        
        oat <- do.call(rbind, lapply(imputed, function(x) {
                rows <- base[c(1, 1), , drop = FALSE]
                rows[[x]] <- rng[x, ]
                enforce_constraints(rows, fixed = c(specified, x))
        }))
        draws <- with_seed(seed, pool[sample.int(nrow(pool), min(N_REF_DRAWS, nrow(pool))), , drop = FALSE])
        joint <- base[rep(1, nrow(draws)), , drop = FALSE]
        joint[, imputed] <- draws[, imputed]
        joint <- enforce_constraints(joint, fixed = specified)
        
        vals  <- tree_values(bm$model, rbind(oat, joint)[, bm$predictors, drop = FALSE])
        med   <- apply(vals, 1, function(z) q7s(sort(z), 0.5))
        n_oat <- nrow(oat)
        
        base_vs    <- sort(v_base)
        base_med   <- q7s(base_vs, 0.5)
        base_range <- q7s(base_vs, RANGE_PROBS)
        eff <- data.table(var = imputed, label = spec_label(imputed),
                          value = as.numeric(unlist(base[1, imputed])), lo = rng[, 1], hi = rng[, 2],
                          median_lo = med[seq(1, n_oat, 2)], median_hi = med[seq(2, n_oat, 2)])
        eff[, `:=`(d_lo = median_lo - base_med, d_hi = median_hi - base_med,
                   impact = abs(median_hi - median_lo))]
        setorder(eff, -impact)
        
        band  <- q7s(sort(med[-seq_len(n_oat)]), RANGE_PROBS)
        ratio <- if (diff(base_range) > 0) diff(band) / diff(base_range) else NA_real_
        mixed <- sort(as.vector(vals[-seq_len(n_oat), , drop = FALSE]))
        
        list(effects = eff, band = band, base_median = base_med, base_range = base_range,
             ratio = ratio, impact = impact_class(ratio), n_draws = nrow(draws), pool_size = nrow(pool),
             mixture = list(values = mixed, median = q7s(mixed, 0.5), range = q7s(mixed, RANGE_PROBS),
                            quantiles = q7s(mixed, taus),
                            pct = if (!is.null(obs)) mid_pct(mixed, obs) else NULL))
}

impact_class <- function(r) {
        if (is.na(r)) return("unknown")
        hit <- names(IMPACT_BREAKS)[r < IMPACT_BREAKS]
        if (length(hit)) hit[1] else "large"
}


# Plots -------------------------------------------------------------------

PAL <- list(fill = "#2e6b8a", line = "#1e7a5a", obs = "#b85c3a", quant = "#4a4a4a",
            bg = "#f5f1ea", border = "#e0dcd4", text = "#4a4a4a", light = "#7a7a7a")

placeholder_plot <- function(msg) {
        ggplot() +
                annotate("text", x = 0.5, y = 0.5, label = msg, size = 4.5, color = PAL$light) +
                theme_void() +
                theme(plot.background = element_rect(fill = "transparent", color = NA))
}

app_plot_theme <- function() {
        theme_minimal(base_family = "sans", base_size = 12) +
                theme(plot.background  = element_rect(fill = "transparent", color = NA),
                      panel.background = element_rect(fill = "transparent", color = NA),
                      panel.grid.major = element_line(color = PAL$border, linewidth = 0.3),
                      panel.grid.minor = element_blank(),
                      axis.text        = element_text(color = PAL$text, size = 9),
                      axis.title       = element_text(color = PAL$text, size = 10),
                      axis.title.y     = element_text(margin = margin(r = 10)),
                      axis.title.x     = element_text(margin = margin(t = 8)),
                      plot.margin      = margin(15, 15, 10, 10))
}

# Legend labels of the two benchmark curves (used in the plot and its explainer)
CURVE_LABELS <- c(fixed  = "Imputed inputs set to single best guesses",
                  varied = "Imputed inputs varied over plausible values")

# qtab: tau, value, lo, hi (whisker = uncertainty interval of the quantile)
build_density_plot <- function(dens, qtab = NULL, observed_value = NULL, mix = NULL,
                               metric_label = "Metric value") {
        if (is.null(dens) || nrow(dens) < 3) return(placeholder_plot("Insufficient data for a density estimate"))
        ymax <- max(c(dens$density, if (!is.null(mix)) mix$density), na.rm = TRUE)
        
        p <- ggplot(dens, aes(x = x, y = density)) +
                geom_area(alpha = 0.18, fill = PAL$fill)
        if (is.null(mix)) {
                p <- p + geom_line(color = PAL$line, linewidth = 0.9)
        } else {
                # Two curves -> legend at the bottom
                curves <- rbind(data.frame(dens[, c("x", "density")], curve = CURVE_LABELS[["fixed"]]),
                                data.frame(mix[, c("x", "density")],  curve = CURVE_LABELS[["varied"]]))
                curves$curve <- factor(curves$curve, levels = unname(CURVE_LABELS))
                p <- p +
                        geom_line(data = curves, aes(x = x, y = density, linetype = curve, color = curve),
                                  linewidth = 0.85) +
                        scale_linetype_manual(values = setNames(c("solid", "dashed"), CURVE_LABELS), name = NULL) +
                        scale_color_manual(values = setNames(c(PAL$line, PAL$quant), CURVE_LABELS), name = NULL)
        }
        
        if (!is.null(qtab) && nrow(qtab)) {
                if (nrow(qtab) >= 2) {
                        shade <- dens[dens$x >= min(qtab$value) & dens$x <= max(qtab$value), ]
                        if (nrow(shade) > 1) p <- p + geom_area(data = shade, alpha = 0.12, fill = PAL$line)
                }
                for (i in seq_len(nrow(qtab))) {
                        med   <- abs(qtab$tau[i] - 0.5) < 0.01
                        col   <- if (med) PAL$line else PAL$quant
                        y_lab <- ymax * if (med) 0.97 else 0.88
                        y_w   <- y_lab - 0.075 * ymax
                        p <- p + geom_vline(xintercept = qtab$value[i], linetype = if (med) "solid" else "dashed",
                                            color = col, linewidth = if (med) 1 else 0.65, alpha = if (med) 0.85 else 0.6)
                        if (all(is.finite(c(qtab$lo[i], qtab$hi[i])))) {
                                tick <- 0.02 * ymax
                                p <- p +
                                        annotate("segment", x = qtab$lo[i], xend = qtab$hi[i], y = y_w, yend = y_w,
                                                 color = col, linewidth = 0.8) +
                                        annotate("segment", x = c(qtab$lo[i], qtab$hi[i]), xend = c(qtab$lo[i], qtab$hi[i]),
                                                 y = y_w - tick, yend = y_w + tick, color = col, linewidth = 0.8)
                        }
                        p <- p + annotate("label", x = qtab$value[i], y = y_lab, label = sprintf("Q%g", qtab$tau[i] * 100),
                                          size = if (med) 3.4 else 2.9, color = col, fill = PAL$bg,
                                          fontface = if (med) "bold" else "plain",
                                          label.padding = unit(0.2, "lines"), label.r = unit(0.2, "lines"),
                                          label.size = 0.2)
                }
        }
        
        if (!is.null(observed_value) && is.finite(observed_value)) {
                p <- p +
                        geom_vline(xintercept = observed_value, color = PAL$obs, linewidth = 1.2) +
                        annotate("label", x = observed_value, y = ymax * 1.08,
                                 label = sprintf("Observed: %s", fmt(observed_value)),
                                 fill = PAL$bg, color = PAL$obs, fontface = "bold", size = 3.2,
                                 label.padding = unit(0.35, "lines"), label.r = unit(0.3, "lines"),
                                 label.size = 0.3)
        }
        
        p + scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
                coord_cartesian(clip = "off") +
                labs(x = metric_label, y = "Density") +
                app_plot_theme() +
                theme(legend.position = "bottom", legend.justification = "left",
                      legend.key.width = unit(2.4, "lines"), legend.text = element_text(color = PAL$text, size = 9),
                      legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(4, "pt"))
}

# Height (px) of the effect plot for n bars
tornado_height <- function(n) 110 + 32 * n

# One bar per imputed input: change in the median benchmark between the lower
# (open point) and upper (filled point) end of its plausible range. Shaded band:
# joint spread over plausible combinations of all imputed inputs.
build_tornado_plot <- function(imp, metric_label, n_show = N_TORNADO) {
        d <- as.data.frame(head(imp$effects, n_show))
        d$display <- sprintf("%s = %s  [%s\u2013%s]", d$label, fmt_in(d$value), fmt_in(d$lo), fmt_in(d$hi))
        d <- d[order(d$impact), ]
        d$display <- factor(d$display, levels = unique(d$display))
        band <- imp$band - imp$base_median
        
        ggplot(d) +
                # Joint spread: shaded band with dashed edges and a label, so it stays
                # recognisable when it is wider than every single bar
                annotate("rect", xmin = band[1], xmax = band[2], ymin = -Inf, ymax = Inf,
                         fill = PAL$line, alpha = 0.10) +
                geom_vline(xintercept = band, color = PAL$line, linewidth = 0.4, linetype = "dashed", alpha = 0.7) +
                annotate("text", x = band[2], y = -Inf, label = "joint spread of all imputed inputs",
                         hjust = 1.03, vjust = -0.6, size = 3, color = PAL$line) +
                geom_vline(xintercept = 0, color = PAL$quant, linewidth = 0.4) +
                geom_segment(aes(x = d_lo, xend = d_hi, y = display, yend = display),
                             linewidth = 4, color = PAL$fill, alpha = 0.45) +
                geom_point(aes(x = d_lo, y = display), shape = 21, fill = "white", color = PAL$fill, size = 2.6, stroke = 0.9) +
                geom_point(aes(x = d_hi, y = display), shape = 21, fill = PAL$fill, color = PAL$fill, size = 2.6) +
                scale_x_continuous(expand = expansion(mult = 0.08)) +
                scale_y_discrete(expand = expansion(add = c(1.1, 0.6))) +   # room for the label below the bars
                labs(x = paste("Change in the median", metric_label), y = NULL) +
                app_plot_theme() +
                theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 9.5))
}



# Figures for files ---------------------------------------------------------

# Formats offered for downloaded figures
FIGURE_FORMATS <- c(PNG = "png", PDF = "pdf")
if (requireNamespace("svglite", quietly = TRUE)) FIGURE_FORMATS <- c(FIGURE_FORMATS, SVG = "svg")

save_plot <- function(p, file, format = "png", width = 9, height = 5.5) {
        if (format == "pdf") {
                ggsave(file, p, width = width, height = height, bg = "white",
                       device = if (capabilities("cairo")) cairo_pdf else "pdf")
        } else {
                ggsave(file, p, width = width, height = height, dpi = 200, bg = "white", device = format,
                       limitsize = FALSE)
        }
}

wrap_text <- function(x, width = 110) paste(strwrap(x, width), collapse = "\n")


# Batch -------------------------------------------------------------------
# Rows are assigned to a model by an optional `metric` column (and an optional
# `taxon` column); without them, the taxon and metric selected in the UI apply.
# Values may be given as display names ("ANOSIM R") or file names ("ANOSIM_R").

choice_key <- function(z) tolower(gsub("[^A-Za-z0-9]+", "", as.character(z)))
match_choice <- function(x, choices, labels = choices) {
        k <- choice_key(x)
        i <- match(k, choice_key(choices))
        j <- match(k, choice_key(labels))
        choices[ifelse(is.na(i), j, i)]
}

find_col <- function(df, name) {
        i <- match(name, tolower(names(df)))
        if (is.na(i)) NULL else names(df)[i]
}

# Labels for rows in figures: first identifier-like column, else the row number
row_labels <- function(df) {
        col <- NULL
        for (nm in c("typology", "id", "name", "label", "site")) {
                col <- find_col(df, nm)
                if (!is.null(col)) break
        }
        lab <- if (is.null(col)) paste("Row", seq_len(nrow(df))) else as.character(df[[col]])
        lab[is.na(lab) | !nzchar(lab)] <- paste("Row", which(is.na(lab) | !nzchar(lab)))
        make.unique(lab, sep = " #")
}

# Model (taxon, metric) per row and a status
resolve_batch_models <- function(df, taxon_default, metric_default) {
        n <- nrow(df)
        tcol <- find_col(df, "taxon"); mcol <- find_col(df, "metric")
        taxa <- unique(MODELS$taxon)
        taxon <- if (is.null(tcol)) rep(taxon_default, n) else
                match_choice(df[[tcol]], taxa, vapply(taxa, function(t) taxon_info(t)$label, ""))
        metric <- rep(NA_character_, n)
        for (t in unique(taxon[!is.na(taxon)])) {
                i  <- which(taxon == t)
                av <- MODELS[taxon == t, metric]
                metric[i] <- if (is.null(mcol)) {
                        if (isTRUE(metric_default %in% av)) metric_default else NA_character_
                } else match_choice(df[[mcol]][i], av, metric_display_label(av))
        }
        status <- ifelse(is.na(taxon), "unknown taxon",
                         ifelse(is.na(metric), if (is.null(mcol)) "no model for the selected metric" else "unknown metric", "ok"))
        list(taxon = taxon, metric = metric, status = status, taxon_col = tcol, metric_col = mcol)
}

# Predictions for the rows of one model
predict_batch_group <- function(bm, dfg, taus, method, k, progress = function(...) NULL) {
        n     <- nrow(dfg)
        preds <- bm$predictors
        cols  <- intersect(preds, names(dfg))
        knn_cache <- new.env(parent = emptyenv())
        rows <- vector("list", n); imputed <- character(n)
        for (i in seq_len(n)) {
                sv  <- as.list(dfg[i, cols, drop = FALSE])
                sv  <- sv[vapply(sv, function(z) is.numeric(z) && is.finite(z), TRUE)]
                inp <- resolve_inputs(bm, sv, method, k, knn_cache)
                rows[[i]]  <- inp$row
                imputed[i] <- paste(inp$imputed, collapse = ";")
        }
        progress(0.5)
        X <- as.data.frame(rbindlist(rows))[, preds, drop = FALSE]
        
        ocol <- find_col(dfg, "observed_value")
        obs  <- if (is.null(ocol)) rep(NA_real_, n) else suppressWarnings(as.numeric(dfg[[ocol]]))
        Q    <- matrix(NA_real_, n, length(taus)); pct <- neff <- rep(NA_real_, n)
        dens <- vector("list", n)
        for (idx in split(seq_len(n), ceiling(seq_len(n) / 500))) {
                tv <- tree_values(bm$model, X[idx, , drop = FALSE])
                for (j in seq_along(idx)) {
                        v  <- tv[j, ]; vs <- sort(v); r <- idx[j]
                        Q[r, ]    <- q7s(vs, taus)
                        dens[[r]] <- density_from_values(v)
                        if (!is.null(bm$training)) neff[r] <- tree_groups(bm, v)$support$n_eff
                        if (is.finite(obs[r])) pct[r] <- mid_pct(vs, obs[r])
                }
        }
        progress(0.5)
        cal <- bm$calibration
        list(imputed = imputed, Q = Q, obs = obs, pct = pct, neff = neff, dens = dens,
             pct_cal = if (is.null(cal)) rep(NA_real_, n) else cal_cdf(cal$cum, cal$edges, pct))
}

run_batch <- function(data, taxon_default, metric_default, taus, method, k,
                      progress = function(amount, detail = NULL) NULL) {
        df <- as.data.frame(data)
        n  <- nrow(df)
        if (!n) stop("The CSV has no rows.")
        res    <- resolve_batch_models(df, taxon_default, metric_default)
        labels <- row_labels(df)
        has_obs <- !is.null(find_col(df, "observed_value"))
        
        imputed <- rep(NA_character_, n); Q <- matrix(NA_real_, n, length(taus))
        pct <- pct_cal <- neff <- rep(NA_real_, n)
        status <- res$status
        keys <- ifelse(status == "ok", paste0(res$taxon, "__", res$metric), NA)
        groups <- list()
        for (key in unique(keys[!is.na(keys)])) {
                idx <- which(keys == key)
                t <- res$taxon[idx[1]]; m <- res$metric[idx[1]]
                progress(0, sprintf("%s, %s: %d rows", taxon_info(t)$label, metric_display_label(m), length(idx)))
                bm <- NULL                                  # release the previous model before loading
                bm <- tryCatch(load_benchmark(t, m), error = function(e) NULL)
                if (is.null(bm)) { status[idx] <- "model could not be loaded"; next }
                g <- predict_batch_group(bm, df[idx, , drop = FALSE], taus, method, k,
                                         function(a) progress(a * length(idx) / n))
                imputed[idx] <- g$imputed; Q[idx, ] <- g$Q
                pct[idx] <- g$pct; pct_cal[idx] <- g$pct_cal; neff[idx] <- g$neff
                groups[[key]] <- list(key = key, taxon = t, metric = m, rows = idx, labels = labels[idx],
                                      taus = taus, Q = g$Q, obs = g$obs, pct = g$pct, pct_cal = g$pct_cal,
                                      dens = g$dens, calibrated = !is.null(bm$calibration))
        }
        
        out <- data.table(df)
        out[, `:=`(model_taxon = ifelse(is.na(res$taxon), NA, vapply(res$taxon, function(t) if (is.na(t)) NA_character_ else taxon_info(t)$label, "")),
                   model_metric = ifelse(is.na(res$metric), NA, metric_display_label(ifelse(is.na(res$metric), "", res$metric))),
                   status = status, imputed_inputs = imputed)]
        for (j in seq_along(taus)) out[, (paste0("Q", taus[j] * 100)) := signif(Q[, j], 4)]
        out[, support_datasets_eff := round(neff, 1)]
        if (has_obs) {
                out[, observed_percentile := round(100 * pct, 1)]
                out[, observed_percentile_calibrated := round(100 * pct_cal, 1)]
        }
        list(table = out, groups = groups, taus = taus, has_obs = has_obs)
}

group_title <- function(g) sprintf("%s - %s", taxon_info(g$taxon)$label, metric_display_label(g$metric))

# Overview of one model: one line per row from the lowest to the highest
# selected quantile, inner line between Q25 and Q75 (if selected), point at the
# median (if selected), observed value as a diamond coloured by its percentile.
build_batch_overview <- function(g, title = TRUE) {
        taus <- g$taus; Q <- g$Q; n <- nrow(Q)
        d <- data.frame(label = factor(g$labels, levels = rev(g$labels)),
                        lo = Q[, which.min(taus)], hi = Q[, which.max(taus)])
        i25 <- match(0.25, taus); i75 <- match(0.75, taus); i50 <- match(0.5, taus)
        p <- ggplot(d, aes(y = label)) +
                geom_segment(aes(x = lo, xend = hi, yend = label), linewidth = 1.2, color = PAL$fill, alpha = 0.45)
        if (!is.na(i25) && !is.na(i75)) {
                p <- p + geom_segment(data = transform(d, a = Q[, i25], b = Q[, i75]),
                                      aes(x = a, xend = b, yend = label), linewidth = 3, color = PAL$fill, alpha = 0.6)
        }
        if (!is.na(i50)) p <- p + geom_point(data = transform(d, m = Q[, i50]), aes(x = m), shape = 21,
                                             fill = "white", color = PAL$line, size = 2.4, stroke = 1)
        ob <- data.frame(label = d$label, obs = g$obs, pct = 100 * g$pct)
        ob <- ob[is.finite(ob$obs), ]
        if (nrow(ob)) {
                ob$class <- cut(ob$pct, c(-Inf, 10, 25, 75, 90, Inf), right = FALSE,
                                labels = c("< 10th", "10th-25th", "25th-75th", "75th-90th", "> 90th"))
                p <- p + geom_point(data = ob, aes(x = obs, fill = class), shape = 23, size = 3, color = "grey20") +
                        scale_fill_manual(values = c("< 10th" = "#b85c3a", "10th-25th" = "#c4a882", "25th-75th" = "#1e7a5a",
                                                     "75th-90th" = "#c4a882", "> 90th" = "#b85c3a"),
                                          drop = FALSE, name = "Observed value,\npercentile")
        }
        sub <- sprintf("Line: Q%g-Q%g%s%s", 100 * min(taus), 100 * max(taus),
                       if (!is.na(i25) && !is.na(i75)) "; thick: Q25-Q75" else "",
                       if (!is.na(i50)) "; circle: median" else "")
        p + labs(x = metric_display_label(g$metric), y = NULL,
                 title = if (title) group_title(g) else NULL, subtitle = sub) +
                app_plot_theme() +
                theme(panel.grid.major.y = element_blank(), legend.position = "right",
                      plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9, color = PAL$light))
}
batch_overview_height <- function(g) min(48, 1.4 + 0.28 * length(g$rows))   # inches

# Predicted distribution of every row of one model, `per_page` panels per page.
# Pages after the first are padded with blank panels, so panels keep their size.
BATCH_PER_PAGE <- 6L
batch_page_height <- function(g) 1.2 + 3.4 * ceiling(min(length(g$rows), BATCH_PER_PAGE) / 2)   # inches
build_batch_density_pages <- function(g, per_page = BATCH_PER_PAGE) {
        idx   <- seq_along(g$rows)
        pages <- split(idx, ceiling(idx / per_page))
        n_first <- length(pages[[1]])
        lapply(pages, function(ii) {
                pad <- if (length(pages) > 1) strrep(" ", seq_len(n_first - length(ii))) else character(0)
                lab <- factor(g$labels[ii], levels = c(g$labels[ii], pad))
                dd  <- rbindlist(lapply(seq_along(ii), function(j) {
                        dj <- g$dens[[ii[j]]]
                        if (is.null(dj)) NULL else data.table(label = lab[j], dj)
                }))
                qd <- rbindlist(lapply(seq_along(ii), function(j)
                        data.table(label = lab[j], tau = g$taus, value = g$Q[ii[j], ])))
                qd[, median := abs(tau - 0.5) < 0.01]
                od <- data.table(label = lab, obs = g$obs[ii], pct = g$pct[ii])[is.finite(obs)]
                p <- ggplot(dd, aes(x = x, y = density)) +
                        geom_area(fill = PAL$fill, alpha = 0.18) +
                        geom_line(color = PAL$line, linewidth = 0.7) +
                        geom_vline(data = qd, aes(xintercept = value, linetype = median), color = PAL$quant,
                                   linewidth = 0.4, show.legend = FALSE) +
                        scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dashed")) +
                        facet_wrap(~ label, ncol = 2, scales = "free", drop = FALSE) +
                        labs(x = metric_display_label(g$metric), y = "Density", title = group_title(g),
                             subtitle = sprintf("Lines: Q%s (solid: median)%s",
                                                paste(100 * g$taus, collapse = ", Q"),
                                                if (nrow(od)) "; red: observed value with its percentile" else "")) +
                        app_plot_theme() +
                        theme(strip.text = element_text(face = "bold", color = PAL$text),
                              plot.title = element_text(face = "bold", size = 12),
                              plot.subtitle = element_text(size = 9, color = PAL$light))
                if (nrow(od)) {
                        p <- p + geom_vline(data = od, aes(xintercept = obs), color = PAL$obs, linewidth = 0.9) +
                                geom_label(data = od, aes(x = obs, y = Inf, label = sprintf("P%.0f", 100 * pct)),
                                           vjust = 1.2, size = 2.8, color = PAL$obs, fill = "white", label.size = 0.2)
                }
                p
        })
}

# Write one CSV, one overview PNG and one distribution PDF per model, plus the combined CSV
write_batch_outputs <- function(batch, dir) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        fwrite(batch$table, file.path(dir, "all_results.csv"))
        pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
        for (g in batch$groups) {
                stem <- g$key
                fwrite(batch$table[g$rows], file.path(dir, paste0(stem, ".csv")))
                
                # Overview: one figure, as PNG and as PDF
                overview <- build_batch_overview(g)
                for (fmt in c("png", "pdf")) {
                        save_plot(overview, file.path(dir, paste0(stem, "_overview.", fmt)), fmt,
                                  width = 8.5, height = batch_overview_height(g))
                }
                
                # Distributions: all pages in one PDF, and one PNG per page
                pages <- build_batch_density_pages(g)
                pdf_dev(file.path(dir, paste0(stem, "_distributions.pdf")), width = 8.27,
                        height = batch_page_height(g), onefile = TRUE)
                for (p in pages) print(p)
                grDevices::dev.off()
                for (i in seq_along(pages)) {
                        png_name <- if (length(pages) == 1) paste0(stem, "_distributions.png")
                        else sprintf("%s_distributions_p%02d.png", stem, i)
                        save_plot(pages[[i]], file.path(dir, png_name), "png",
                                  width = 8.27, height = batch_page_height(g))
                }
        }
        invisible(list.files(dir))
}

# zip if a zip program is available, otherwise tar.gz (pure R)
ARCHIVE_EXT <- if (nzchar(Sys.which("zip"))) "zip" else "tar.gz"
make_archive <- function(dir, file) {
        if (ARCHIVE_EXT == "zip") {
                utils::zip(file, list.files(dir, full.names = TRUE), flags = "-q9Xj")   # -j: no folder paths
        } else {
                old <- setwd(dir); on.exit(setwd(old))
                utils::tar(file, list.files("."), compression = "gzip", tar = "internal")
        }
}


# Theme -------------------------------------------------------------------
custom_theme <- bs_theme(
        version = 5,
        bg = "#f5f1ea",
        fg = "#1a1a1a",
        primary = "#1e7a5a",
        secondary = "#4a4a4a",
        success = "#2a9d6e",
        info = "#2e6b8a",
        warning = "#c4a882",
        danger = "#b85c3a",
        base_font = font_google("Source Sans 3"),
        heading_font = font_google("Instrument Serif"),
        code_font = font_google("JetBrains Mono"),
        font_scale = 0.95,
        `enable-rounded` = TRUE
) |>
        bs_add_rules("
    :root {
      --river-deep: #0a1628;
      --river-dark: #122240;
      --water-blue: #2e6b8a;
      --water-light: #4a9aba;
      --foam: #b8d8e8;
      --accent: #1e7a5a;
      --accent-light: #2a9d6e;
      --parchment: #f5f1ea;
      --white: #fefefe;
      --text-primary: #1a1a1a;
      --text-secondary: #4a4a4a;
      --text-light: #7a7a7a;
      --border: #e0dcd4;
      --card-shadow: 0 2px 12px rgba(10, 22, 40, 0.06);
    }

    body {
      background: var(--parchment) !important;
      -webkit-font-smoothing: antialiased;
    }

    /* ─── NAVBAR ─── */
    .navbar {
      background: rgba(245, 241, 234, 0.95) !important;
      backdrop-filter: blur(12px);
      -webkit-backdrop-filter: blur(12px);
      border-bottom: 1px solid var(--border) !important;
      box-shadow: none !important;
    }

    .navbar-brand {
      font-family: 'Instrument Serif', Georgia, serif !important;
      font-size: 1.3rem !important;
      color: var(--river-deep) !important;
      letter-spacing: -0.02em;
    }

    .navbar .nav-link {
      color: var(--text-secondary) !important;
      font-size: 0.85rem;
      font-weight: 400;
      letter-spacing: 0.03em;
      text-transform: uppercase;
      transition: color 0.2s ease;
    }

    .navbar .nav-link:hover,
    .navbar .nav-link.active {
      color: var(--accent) !important;
    }

    .navbar .nav-link.active {
      border-bottom: 2px solid var(--accent);
    }

    /* ─── CARDS ─── */
    .card {
      background: var(--white);
      border: 1px solid var(--border);
      border-radius: 12px;
      box-shadow: var(--card-shadow);
      transition: transform 0.2s ease, box-shadow 0.2s ease;
    }

    .card:hover {
      transform: translateY(-1px);
      box-shadow: 0 4px 20px rgba(10, 22, 40, 0.08);
    }

    .card-header {
      background: transparent !important;
      border-bottom: 1px solid var(--border);
      font-family: 'Instrument Serif', Georgia, serif;
      font-weight: 400;
      font-size: 1.1rem;
      color: var(--river-deep);
      letter-spacing: -0.01em;
      padding: 1rem 1.25rem;
    }

    /* ─── SIDEBAR ─── */
    .sidebar {
      background: var(--white) !important;
      border-right: 1px solid var(--border) !important;
    }

    .sidebar .sidebar-title {
      font-family: 'Instrument Serif', Georgia, serif !important;
      color: var(--river-deep) !important;
      font-size: 1.2rem;
    }

    /* ─── FORM CONTROLS ─── */
    .form-control, .form-select {
      background: var(--parchment) !important;
      border: 1px solid var(--border) !important;
      color: var(--text-primary) !important;
      border-radius: 8px;
      transition: all 0.2s ease;
      font-size: 0.9rem;
    }

    .form-control:focus, .form-select:focus {
      background: var(--white) !important;
      border-color: var(--accent) !important;
      box-shadow: 0 0 0 3px rgba(30, 122, 90, 0.12) !important;
    }

    .form-label, .control-label {
      font-size: 0.85rem;
      font-weight: 500;
      color: var(--text-secondary);
      letter-spacing: 0.01em;
    }

    /* ─── BUTTONS ─── */
    .btn-primary {
      background: linear-gradient(135deg, var(--accent), var(--water-blue)) !important;
      border: none !important;
      font-weight: 500;
      letter-spacing: 0.02em;
      padding: 0.7rem 1.8rem;
      border-radius: 8px;
      transition: all 0.3s ease;
    }

    .btn-primary:hover {
      transform: translateY(-1px);
      box-shadow: 0 4px 16px rgba(30, 122, 90, 0.25);
    }

    .btn-outline-secondary {
      border-color: var(--border) !important;
      color: var(--text-secondary) !important;
      border-radius: 8px;
    }

    .btn-outline-secondary:hover {
      background: var(--parchment) !important;
      border-color: var(--water-blue) !important;
      color: var(--water-blue) !important;
    }

    /* ─── NAV PILLS (parameter tabs) ─── */
    .nav-pills .nav-link {
      color: var(--text-secondary);
      border-radius: 8px;
      margin: 2px;
      font-size: 0.82rem;
      font-weight: 500;
      transition: all 0.2s ease;
    }

    .nav-pills .nav-link:hover {
      background: var(--parchment);
      color: var(--accent);
    }

    .nav-pills .nav-link.active {
      background: linear-gradient(135deg, var(--accent), var(--water-blue)) !important;
      color: white !important;
    }

    /* ─── PREDICTION DISPLAY ─── */
    .prediction-value {
      font-family: 'Instrument Serif', Georgia, serif;
      font-size: 3.5rem;
      font-weight: 400;
      background: linear-gradient(135deg, var(--accent), var(--water-blue));
      -webkit-background-clip: text;
      -webkit-text-fill-color: transparent;
      background-clip: text;
      letter-spacing: -0.03em;
    }

    .quantile-badge {
      display: inline-block;
      padding: 0.45rem 0.9rem;
      border-radius: 20px;
      font-weight: 500;
      font-size: 0.85rem;
      margin: 0.2rem;
    }

    .quantile-low {
      background: rgba(184, 92, 58, 0.1);
      color: #b85c3a;
    }

    .quantile-mid {
      background: rgba(46, 107, 138, 0.1);
      color: var(--water-blue);
    }

    .quantile-high {
      background: rgba(30, 122, 90, 0.1);
      color: var(--accent);
    }

    /* ─── PARAMETER GROUP HEADERS ─── */
    .parameter-group-title {
      font-family: 'Source Sans 3', sans-serif;
      font-size: 0.72rem;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.12em;
      color: var(--accent);
      margin-bottom: 0.8rem;
      padding-bottom: 0.5rem;
      border-bottom: 1px solid var(--border);
    }

    /* ─── CHECKBOXES & RADIOS ─── */
    .form-check-input:checked {
      background-color: var(--accent) !important;
      border-color: var(--accent) !important;
    }

    .form-check-input:focus {
      box-shadow: 0 0 0 3px rgba(30, 122, 90, 0.15) !important;
    }

    .form-check-label {
      font-size: 0.88rem;
      color: var(--text-secondary);
    }

    /* ─── STATUS INDICATOR ─── */
    .status-indicator {
      width: 8px;
      height: 8px;
      border-radius: 50%;
      display: inline-block;
      margin-right: 8px;
      animation: pulse 2s infinite;
    }

    .status-ready {
      background: var(--accent);
    }

    .status-loading {
      background: #c4a882;
    }

    @keyframes pulse {
      0%, 100% { opacity: 1; }
      50% { opacity: 0.5; }
    }

    /* ─── METRIC / INFO CARDS ─── */
    .metric-card {
      text-align: center;
      padding: 1.25rem;
    }

    .metric-label {
      font-size: 0.72rem;
      text-transform: uppercase;
      letter-spacing: 0.1em;
      color: var(--text-light);
      margin-bottom: 0.4rem;
      font-weight: 500;
    }

    /* ─── BADGES ─── */
    .badge.bg-success {
      background: rgba(30, 122, 90, 0.12) !important;
      color: var(--accent) !important;
      font-weight: 500;
    }

    .badge.bg-warning {
      background: rgba(196, 168, 130, 0.2) !important;
      color: #8a6d3b !important;
      font-weight: 500;
    }

    .badge.bg-info {
      background: rgba(46, 107, 138, 0.12) !important;
      color: var(--water-blue) !important;
      font-weight: 500;
    }

    /* ─── DATA TABLES ─── */
    .dataTables_wrapper {
      color: var(--text-primary) !important;
    }

    .dataTables_wrapper .dataTables_filter input,
    .dataTables_wrapper .dataTables_length select {
      background: var(--parchment) !important;
      border: 1px solid var(--border) !important;
      color: var(--text-primary) !important;
      border-radius: 6px;
    }

    table.dataTable {
      border-collapse: collapse !important;
    }

    table.dataTable thead th {
      background: var(--parchment) !important;
      border-bottom: 2px solid var(--accent) !important;
      color: var(--river-deep);
      font-family: 'Source Sans 3', sans-serif;
      font-weight: 600;
      font-size: 0.85rem;
      letter-spacing: 0.02em;
    }

    table.dataTable tbody tr {
      background: transparent !important;
      border-bottom: 1px solid var(--border);
    }

    table.dataTable tbody tr:hover {
      background: rgba(46, 107, 138, 0.04) !important;
    }

    table.dataTable tbody td {
      color: var(--text-secondary);
      font-size: 0.88rem;
    }

    /* ─── MARKDOWN CONTENT (About page) ─── */
    .card h3 {
      font-family: 'Instrument Serif', Georgia, serif;
      font-weight: 400;
      color: var(--river-deep);
      font-size: 1.25rem;
      margin-top: 1.5rem;
    }

    .card code {
      font-size: 0.85em;
      color: var(--water-blue);
      background: rgba(46, 107, 138, 0.06);
      padding: 0.15em 0.4em;
      border-radius: 4px;
    }

    .card pre {
      background: var(--parchment) !important;
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 1rem;
    }

    .card pre code {
      background: transparent;
      padding: 0;
      font-size: 0.82rem;
    }

    .card ul li {
      color: var(--text-secondary);
      margin-bottom: 0.3rem;
    }

    /* ─── SCROLLBARS ─── */
    ::-webkit-scrollbar {
      width: 6px;
    }

    ::-webkit-scrollbar-track {
      background: var(--parchment);
    }

    ::-webkit-scrollbar-thumb {
      background: var(--border);
      border-radius: 3px;
    }

    ::-webkit-scrollbar-thumb:hover {
      background: var(--water-light);
    }

    /* ─── NOTIFICATION TOASTS ─── */
    .shiny-notification {
      background: var(--white);
      border: 1px solid var(--border);
      color: var(--text-primary);
      border-radius: 10px;
      box-shadow: 0 4px 20px rgba(10, 22, 40, 0.1);
    }

    /* ─── HORIZONTAL RULES ─── */
    hr {
      border-color: var(--border);
      opacity: 0.6;
    }

    /* ─── FILE INPUT ─── */
    .form-control[type='file'] {
      font-size: 0.85rem;
    }

    /* ─── IMPUTED VALUES HIGHLIGHT ─── */
    .imputed-summary {
      background: rgba(46, 107, 138, 0.04);
      border: 1px solid rgba(46, 107, 138, 0.15);
      border-radius: 8px;
      padding: 0.75rem 1rem;
      margin-top: 0.75rem;
      font-size: 0.82rem;
    }

    .imputed-summary .imputed-header {
      font-weight: 600;
      color: var(--water-blue);
      font-size: 0.72rem;
      text-transform: uppercase;
      letter-spacing: 0.08em;
      margin-bottom: 0.4rem;
    }

    .imputed-var {
      display: inline-block;
      background: rgba(46, 107, 138, 0.08);
      border-radius: 4px;
      padding: 0.15rem 0.5rem;
      margin: 0.15rem;
      font-family: 'JetBrains Mono', monospace;
      font-size: 0.78rem;
      color: var(--text-secondary);
    }

    .imputed-var .imputed-range {
      color: var(--text-light);
    }

    /* ─── DENSITY PLOT CONTAINER ─── */
    .density-plot-container {
      border-top: 1px solid var(--border);
      margin-top: 1.5rem;
      padding-top: 1rem;
    }

    .density-plot-header {
      font-family: 'Instrument Serif', Georgia, serif;
      font-size: 1rem;
      color: var(--river-deep);
      margin-bottom: 0.5rem;
    }

    .density-plot-subtitle {
      font-size: 0.78rem;
      color: var(--text-light);
      margin-bottom: 0.75rem;
    }

    /* ─── OBSERVED VALUE INPUT STYLING ─── */
    .observed-value-section {
      background: rgba(184, 92, 58, 0.04);
      border: 1px solid rgba(184, 92, 58, 0.15);
      border-radius: 8px;
      padding: 0.75rem 1rem;
      margin-top: 0.5rem;
    }

    .observed-value-section .control-label {
      color: #b85c3a !important;
      font-weight: 600;
    }

    /* ─── PERCENTILE BADGE (for observed value) ─── */
    .percentile-badge {
      display: inline-block;
      padding: 0.3rem 0.7rem;
      border-radius: 16px;
      font-weight: 600;
      font-size: 0.82rem;
      background: rgba(184, 92, 58, 0.1);
      color: #b85c3a;
      margin-top: 0.3rem;
    }

    /* ─── TWO-CURVE EXPLAINER (below the density plot) ─── */
    .curve-explainer {
      background: rgba(46, 107, 138, 0.04);
      border: 1px solid rgba(46, 107, 138, 0.15);
      border-radius: 8px;
      padding: 0.75rem 1rem;
      margin-top: 0.5rem;
      font-size: 0.84rem;
      color: var(--text-secondary);
    }

    .curve-explainer .curve-row {
      display: flex;
      gap: 0.6rem;
      align-items: baseline;
      margin-bottom: 0.35rem;
    }

    .curve-key {
      flex: 0 0 28px;
      height: 0;
      border-top: 3px solid #1e7a5a;
      transform: translateY(-0.25em);
    }

    .curve-key.dashed {
      border-top: 2px dashed #4a4a4a;
    }

    .curve-explainer .curve-reading {
      margin-top: 0.5rem;
      padding-top: 0.5rem;
      border-top: 1px solid rgba(46, 107, 138, 0.15);
      color: var(--text-light);
    }

    /* ─── UNCERTAINTY TABLE ─── */
    .uncert-table {
      margin: 0.75rem auto 0;
      border-collapse: collapse;
      font-size: 0.85rem;
      min-width: 60%;
    }

    .uncert-table th {
      font-weight: 600;
      font-size: 0.78rem;
      color: var(--river-deep);
      border-bottom: 2px solid var(--accent);
      padding: 0.35rem 0.75rem;
      text-align: right;
    }

    .uncert-table td {
      padding: 0.3rem 0.75rem;
      border-bottom: 1px solid var(--border);
      color: var(--text-secondary);
      text-align: right;
      font-variant-numeric: tabular-nums;
    }

    .uncert-table th:first-child, .uncert-table td:first-child {
      text-align: left;
    }

    .uncert-ci {
      color: var(--text-light);
      font-size: 0.78rem;
    }

    .uncert-note {
      font-size: 0.75rem;
      color: var(--text-light);
      max-width: 680px;
      margin: 0.5rem auto 0;
      text-align: left;
    }

    /* ─── IMPUTATION IMPACT ─── */
    .impact-badge {
      display: inline-block;
      padding: 0.25rem 0.7rem;
      border-radius: 14px;
      font-weight: 600;
      font-size: 0.78rem;
      text-transform: uppercase;
      letter-spacing: 0.06em;
      margin-bottom: 0.5rem;
    }

    .impact-negligible { background: rgba(30, 122, 90, 0.12); color: var(--accent); }
    .impact-minor      { background: rgba(46, 107, 138, 0.12); color: var(--water-blue); }
    .impact-moderate   { background: rgba(196, 168, 130, 0.25); color: #8a6d3b; }
    .impact-large      { background: rgba(184, 92, 58, 0.12); color: #b85c3a; }
    .impact-unknown    { background: var(--parchment); color: var(--text-light); }

    .imputation-text p {
      font-size: 0.9rem;
      color: var(--text-secondary);
      margin-bottom: 0.4rem;
    }
  ")


# UI ----------------------------------------------------------------------

quantile_choices <- c("5%" = 0.05, "10%" = 0.10, "25%" = 0.25, "50% (Median)" = 0.50,
                      "75%" = 0.75, "90%" = 0.90, "95%" = 0.95)
impute_choices <- c("Conditional averages (KNN)" = "conditional", "Training medians" = "marginal")
taxa           <- unique(MODELS$taxon)
taxon_choices  <- if (length(taxa)) setNames(taxa, vapply(taxa, function(t) taxon_info(t)$label, "")) else c("No models found" = "")
taxon_default  <- if ("macroinvertebrates" %in% taxa) "macroinvertebrates" else unname(taxon_choices[1])
no_models_msg  <- if (!length(taxa)) div(class = "small mb-2", style = "color: var(--danger);", icon("triangle-exclamation"),
                                         sprintf(" No model files found in %s", MODEL_DIR))

ui <- page_navbar(
        title = tags$span(tags$span(class = "status-indicator status-ready"), "QRF Predictor"),
        theme = custom_theme,
        fillable = TRUE,
        
        # Single prediction ------------------------------------------------
        nav_panel(
                title = "Predict",
                icon = icon("chart-line"),
                
                layout_sidebar(
                        sidebar = sidebar(
                                title = "Parameters",
                                width = 380,
                                
                                card(
                                        card_header("Model Selection"),
                                        no_models_msg,
                                        radioButtons("taxon_group", "Taxonomic Group", choices = taxon_choices,
                                                     selected = taxon_default),
                                        hr(),
                                        selectInput("metric", "Evaluation Metric", choices = NULL, selected = NULL),
                                        checkboxGroupInput("quantiles", "Prediction Quantiles", choices = quantile_choices,
                                                           selected = c(0.10, 0.50, 0.90), inline = TRUE),
                                        hr(),
                                        div(
                                                class = "observed-value-section",
                                                numericInput("observed_value",
                                                             tags$span(icon("crosshairs"), " Your Observed Value"),
                                                             value = NA, step = 0.001),
                                                tags$small(class = "text-muted",
                                                           "Optional. Enter your typology's actual metric value to see ",
                                                           "where it falls in the predicted distribution.")
                                        ),
                                        hr(),
                                        radioButtons("impute_method", "Unspecified Parameters",
                                                     choices = impute_choices, selected = "conditional"),
                                        conditionalPanel(
                                                condition = "input.impute_method == 'conditional'",
                                                sliderInput("knn_k", "Neighbors (k)", min = 5, max = 200, value = 50, step = 5)
                                        ),
                                        uiOutput("reference_status")
                                ),
                                
                                card(
                                        card_header("Active Parameters"),
                                        p(class = "text-muted small",
                                          "Toggle which parameters you specify. Unselected parameters are imputed ",
                                          "with the method chosen above."),
                                        div(class = "small mb-2",
                                            actionLink("toggle_all", "Select all"), " \u00b7 ",
                                            actionLink("toggle_none", "Select none")),
                                        do.call(navset_pill, c(list(id = "param_tabs"), lapply(PREDICTOR_GROUPS, function(g)
                                                nav_panel(g, div(class = "mt-3", uiOutput(paste0("toggles_", g)))))))
                                ),
                                
                                div(
                                        class = "d-grid gap-2 mt-3",
                                        actionButton("predict", "Generate Prediction", class = "btn-primary btn-lg",
                                                     icon = icon("wand-magic-sparkles")),
                                        actionButton("reset", "Reset to Training Medians", class = "btn-outline-secondary",
                                                     icon = icon("rotate-left"))
                                )
                        ),
                        
                        # fill = FALSE: cards keep their natural height and the main area scrolls.
                        # In a filled layout the cards share the window height, which squeezes
                        # or clips the plots.
                        layout_columns(
                                col_widths = c(12, 12, 8, 4),
                                fill = FALSE,
                                
                                card(
                                        card_header(class = "d-flex justify-content-between align-items-center",
                                                    "Prediction Result", uiOutput("model_status")),
                                        uiOutput("prediction_display"),
                                        div(
                                                class = "density-plot-container",
                                                div(class = "density-plot-header", icon("chart-area"), " Predicted Distribution"),
                                                uiOutput("density_subtitle"),
                                                plotOutput("density_plot", height = "330px", fill = FALSE),
                                                div(class = "d-flex justify-content-end align-items-center gap-3 mt-1",
                                                    radioButtons("figure_format", NULL, choices = FIGURE_FORMATS,
                                                                 selected = "png", inline = TRUE),
                                                    downloadButton("download_density", "Download figure",
                                                                   class = "btn-outline-secondary btn-sm")),
                                                uiOutput("observed_percentile_display"),
                                                uiOutput("curve_explainer")
                                        )
                                ),
                                
                                card(
                                        card_header(icon("sliders"), " Imputed Inputs and Their Effect"),
                                        uiOutput("imputation_panel")
                                ),
                                
                                card(
                                        card_header("Parameter Values"),
                                        div(style = "max-height: 500px; overflow-y: auto;", uiOutput("parameter_inputs"))
                                ),
                                
                                card(card_header("Model Info"), uiOutput("model_info"))
                        )
                )
        ),
        
        # Batch prediction -------------------------------------------------
        nav_panel(
                title = "Batch",
                icon = icon("table"),
                
                layout_columns(
                        col_widths = c(4, 8),
                        fill = FALSE,
                        
                        card(
                                card_header("Upload Data"),
                                fileInput("batch_file", "CSV File", accept = c(".csv", ".CSV")),
                                tags$small(class = "text-muted",
                                           "One row per typology; predictor columns named as in the Parameter Reference ",
                                           "(About tab). Missing columns or empty cells are imputed. Optional columns: ",
                                           tags$code("metric"), " (model per row, so one file can hold several metrics), ",
                                           tags$code("taxon"), ", ", tags$code("observed_value"), " (adds percentiles) and ",
                                           tags$code("typology"), " (row names in the figures)."),
                                radioButtons("batch_taxon_group", "Taxonomic Group", choices = taxon_choices,
                                             selected = taxon_default, inline = TRUE),
                                tags$small(class = "text-muted d-block mb-2", "Used for rows without a ", tags$code("taxon"), " column."),
                                selectInput("batch_metric", "Evaluation Metric", choices = NULL, selected = NULL),
                                tags$small(class = "text-muted d-block mb-2", "Used only if the file has no ", tags$code("metric"), " column."),
                                checkboxGroupInput("batch_quantiles", "Quantiles to Predict",
                                                   choices = c("10%" = 0.1, "25%" = 0.25, "50%" = 0.5, "75%" = 0.75, "90%" = 0.9),
                                                   selected = c(0.1, 0.25, 0.5, 0.75, 0.9)),
                                hr(),
                                radioButtons("batch_impute_method", "Missing Value Handling",
                                             choices = impute_choices, selected = "conditional"),
                                conditionalPanel(
                                        condition = "input.batch_impute_method == 'conditional'",
                                        sliderInput("batch_knn_k", "Neighbors (k)", min = 5, max = 200, value = 50, step = 5)
                                ),
                                uiOutput("batch_column_status"),
                                div(class = "d-grid gap-2 mt-3",
                                    actionButton("batch_predict", "Run Batch Prediction", class = "btn-primary", icon = icon("play")),
                                    downloadButton("download_batch_archive",
                                                   sprintf("Download results and figures (.%s)", ARCHIVE_EXT),
                                                   class = "btn-outline-secondary"),
                                    downloadButton("download_results", "Download combined CSV", class = "btn-outline-secondary")),
                                tags$small(class = "text-muted mt-2",
                                           "The archive holds, per metric, a CSV, an overview figure (PNG) and the predicted ",
                                           "distribution of every row (PDF), plus all_results.csv.")
                        ),
                        
                        div(
                                card(card_header("Results"), DTOutput("batch_results", fill = FALSE)),
                                card(card_header(icon("chart-simple"), " Figures"),
                                     uiOutput("batch_figure_panel"))
                        )
                )
        ),
        
        # About ------------------------------------------------------------
        nav_panel(
                title = "About",
                icon = icon("circle-info"),
                
                layout_columns(
                        col_widths = c(6, 6),
                        
                        card(
                                card_header("About This App"),
                                markdown("
This application provides predictions from **Quantile Random Forest** models
fitted on freshwater typology evaluation metrics: Classification Strength,
ANOSIM R, AucZeta, PERMANOVA R2, PERMANOVA Fuzzy R2 and Fuzzy Mantel. For a typology described by
its inputs, the model returns the distribution of the metric expected from
simulated typologies with the same characteristics: the benchmark.

### Taxonomic Groups

- 🌿 **Diatoms** — Microscopic algae
- 🐟 **Fishes** — Fish assemblages
- 🦗 **Macroinvertebrates** — Benthic invertebrates
- 🌱 **Macrophytes** — Aquatic plants

### How It Works

1. Select a taxonomic group and an evaluation metric
2. Toggle the inputs you can specify and enter their values
3. Choose how the other inputs are imputed:
   - *Conditional averages (KNN)*: averages over the *k* training typologies
     closest to your specified inputs (standardised Euclidean distance)
   - *Training medians*: medians over all training typologies
4. Optionally enter your observed metric value
5. Click predict

### The Predicted Distribution

Each of the forest's trees contributes one prediction. Their distribution is
the benchmark; the plot shows its kernel density, the selected quantiles and,
if given, your observed value with its percentile.

### Uncertainty of Quantiles and Percentiles

Intervals are 90% intervals.

- **Estimate**: each tree predicts one value from the training data, and the
  app traces that value to the training data set it came from. Resampling the
  trees and these data sets shows how much a quantile, or the percentile of
  your observed value, depends on which data sets were available for training.
  The forest structure is held fixed, so the interval does not include the
  variation a refit would add.
- **Data support**: the number of training data sets the trees' predictions
  come from, and their effective number. Few data sets mean a fragile
  benchmark.
- **Calibrated** (when the model's held-out results are available): quantiles
  are matched to typologies from data sets held out during fitting (test set
  and out-of-fold cross-validation). The calibrated Q10, for example, is the
  value below which 10% of held-out typologies fell. Calibration is global: it
  corrects the average reliability of the quantiles, not their reliability for
  one combination of inputs.

Inputs outside the training range are flagged: the forest does not
extrapolate, so predictions stay flat beyond the range. None of the intervals
covers uncertainty in your own input values.

### Imputed Inputs

Every imputed input gets a plausible range: the 10th–90th percentile among the
training typologies used for imputation (the *k* nearest neighbours, or all
training typologies). The app then

- varies each imputed input over its range, one at a time, and shows the change
  in the median benchmark (bar chart);
- draws up to 200 of these typologies and uses their values for all imputed
  inputs together, which keeps their covariance. The spread of the median over
  these draws, relative to the width of the benchmark's own 80% interval,
  classifies the impact of imputation (negligible < 10% ≤ minor < 33% ≤
  moderate < 66% ≤ large);
- pools the predictions over these draws into a benchmark that accounts for
  imputation (dashed line in the plot).

### Files

Models are read from the folder `models` next to the app, or the folder set in
the environment variable `MODEL_DIR`: `<taxon>__<metric>.rds` (the forest and its
prepped recipe, which contains the training typologies) and, where present,
`<taxon>__<metric>_intervals.rds` (held-out evaluation, used for calibration).

### Online and Full Models

The online app serves reduced forests: the first 1,000 of each fitted forest's
3,000 trees, unchanged. Model Info shows the number of trees a model keeps.
Fewer trees add Monte Carlo noise to quantiles and percentiles, which the
intervals include. The full forests can be downloaded and run in this app
locally; they need more memory than the reduced ones.
        ")
                        ),
                        
                        card(
                                card_header("Parameter Reference"),
                                div(style = "max-height: 600px; overflow-y: auto;", DTOutput("param_reference"))
                        )
                )
        )
)


# Server ------------------------------------------------------------------

server <- function(input, output, session) {
        
        # The session keeps the model's metadata only (benchmark_light); the forest stays
        # in the shared cache and is fetched when a prediction runs, reloaded if it was
        # released to make room for another model.
        benchmark         <- reactiveVal(NULL)   # benchmark_light(load_benchmark())
        benchmark_full    <- function() {
                meta <- benchmark()
                if (is.null(meta)) return(NULL)
                if (!is.null(.memo$items[[meta$key]])) return(load_benchmark(meta$taxon, meta$metric))
                withProgress(message = "Reloading model", detail = "It was released to make room for another model...",
                             value = 0.3, load_benchmark(meta$taxon, meta$metric))
        }
        prediction_result <- reactiveVal(NULL)
        batch_results     <- reactiveVal(NULL)
        
        # ---- Metric lists ----------------------------------------------------
        update_metrics <- function(input_id, taxon) {
                metrics <- get_available_metrics(taxon)
                current <- isolate(input[[input_id]])
                updateSelectInput(session, input_id, choices = metrics,
                                  selected = if (isTRUE(current %in% metrics)) current else unname(metrics[1]))
        }
        observeEvent(input$taxon_group,       update_metrics("metric", input$taxon_group))
        observeEvent(input$batch_taxon_group, update_metrics("batch_metric", input$batch_taxon_group))
        
        # ---- Load model ------------------------------------------------------
        observeEvent(list(input$taxon_group, input$metric), {
                req(input$taxon_group, input$metric, nzchar(input$metric))
                bm <- withProgress(message = "Loading model", detail = "Reading model file (first use can take a while)...",
                                   value = 0.3, {
                                           tryCatch(load_benchmark(input$taxon_group, input$metric), error = function(e) {
                                                   showNotification(conditionMessage(e), type = "error", duration = 8)
                                                   NULL
                                           })
                                   })
                benchmark(benchmark_light(bm))
                prediction_result(NULL)
                if (!is.null(bm)) {
                        showNotification(paste("Model loaded:", metric_display_label(input$metric)), type = "message", duration = 3)
                        if (is.null(bm$training)) {
                                showNotification("No training data in this model file: unspecified inputs use hard-coded defaults.",
                                                 type = "warning", duration = 6)
                        }
                }
        })
        
        # ---- Status displays -------------------------------------------------
        output$reference_status <- renderUI({
                bm <- benchmark()
                if (is.null(bm)) return(NULL)
                tr <- bm$training
                if (!is.null(tr)) {
                        div(class = "small mt-2", style = "color: var(--accent);", icon("check-circle"),
                            sprintf(" Training typologies loaded (%s from %d data sets)",
                                    format(tr$n_rows, big.mark = ","), tr$n_datasets))
                } else {
                        div(class = "small mt-2", style = "color: var(--danger);", icon("triangle-exclamation"),
                            " The model file carries no training data: unspecified inputs use hard-coded defaults, ",
                            "without ranges, imputation effects or data-set support.")
                }
        })
        
        output$model_status <- renderUI({
                bm <- benchmark()
                if (!is.null(bm)) {
                        tags$span(class = "badge bg-success", icon("check"),
                                  paste(taxon_info(bm$taxon)$label, "- Ready"))
                } else {
                        tags$span(class = "badge bg-warning", icon("exclamation-triangle"), "No Model")
                }
        })
        
        # ---- Toggles and inputs ----------------------------------------------
        params_in_group <- function(group) {
                p  <- names(PREDICTOR_SPEC)[vapply(PREDICTOR_SPEC, function(x) x$group == group, TRUE)]
                bm <- benchmark()
                if (is.null(bm)) p else intersect(p, bm$predictors)
        }
        
        param_label <- function(p) {
                spec <- PREDICTOR_SPEC[[p]]
                tooltip(tags$span(spec$label, " ", icon("circle-info", style = "color: var(--text-light);")),
                        spec$help, placement = "right")
        }
        
        lapply(PREDICTOR_GROUPS, function(g) {
                output[[paste0("toggles_", g)]] <- renderUI({
                        params <- params_in_group(g)
                        if (!length(params)) return(p(class = "text-muted small", "Not used by this model."))
                        tagList(lapply(params, function(p) {
                                id <- paste0("toggle_", p)
                                checkboxInput(id, PREDICTOR_SPEC[[p]]$label, value = isolate(input[[id]]) %||% TRUE)
                        }))
                })
                # Render toggles on hidden tabs too, so "Select all/none" reaches every parameter
                outputOptions(output, paste0("toggles_", g), suspendWhenHidden = FALSE)
        })
        
        set_toggles <- function(value) {
                for (p in names(PREDICTOR_SPEC)) updateCheckboxInput(session, paste0("toggle_", p), value = value)
        }
        observeEvent(input$toggle_all,  set_toggles(TRUE))
        observeEvent(input$toggle_none, set_toggles(FALSE))
        
        is_toggled <- function(p) isTRUE(input[[paste0("toggle_", p)]] %||% TRUE)
        
        output$parameter_inputs <- renderUI({
                groups <- lapply(PREDICTOR_GROUPS, function(group) {
                        params <- Filter(is_toggled, params_in_group(group))
                        if (!length(params)) return(NULL)
                        div(
                                div(class = "parameter-group-title", group),
                                lapply(params, function(p) {
                                        spec <- PREDICTOR_SPEC[[p]]
                                        id   <- paste0("param_", p)
                                        numericInput(id, label = param_label(p),
                                                     value = isolate(input[[id]]) %||% spec$value,
                                                     min = spec$min, max = spec$max, step = spec$step %||% 1)
                                })
                        )
                })
                groups <- Filter(Negate(is.null), groups)
                if (!length(groups)) {
                        return(p(class = "text-muted", "No parameters selected: all inputs will be imputed."))
                }
                tagList(groups)
        })
        
        observeEvent(input$reset, {
                med <- benchmark()$training$medians
                for (p in names(PREDICTOR_SPEC)) {
                        val <- if (!is.null(med) && p %in% names(med)) med[[p]] else PREDICTOR_SPEC[[p]]$value
                        if (is_int_var(p)) val <- round(val)
                        updateNumericInput(session, paste0("param_", p), value = signif(val, 4))
                }
                updateNumericInput(session, "observed_value", value = NA)
                prediction_result(NULL)
                showNotification(if (is.null(med)) "Parameters reset to defaults" else "Parameters reset to training medians",
                                 type = "message")
        })
        
        # ---- Model info ------------------------------------------------------
        output$model_info <- renderUI({
                bm <- benchmark()
                if (is.null(bm)) {
                        return(div(class = "text-center p-4", style = "color: var(--text-light);",
                                   icon("folder-open", class = "fa-3x mb-3"), p("No model loaded"),
                                   p(class = "small", "Select a taxonomic group and metric")))
                }
                ti <- taxon_info(bm$taxon)
                tr    <- bm$training
                cov90 <- if (!is.null(bm$coverage)) bm$coverage[abs(bm$coverage$nominal - 0.9) < 1e-9] else NULL
                info_row <- function(label, value) div(class = "metric-card",
                                                       div(class = "metric-label", label),
                                                       div(class = "h4", style = "color: var(--river-deep); font-family: 'Instrument Serif', serif;", value))
                div(
                        div(class = "text-center mb-3 p-2", style = "background: rgba(46, 107, 138, 0.06); border-radius: 8px;",
                            icon(ti$icon, class = "fa-2x mb-2", style = "color: var(--water-blue);"),
                            div(class = "fw-bold", style = "color: var(--river-deep);", ti$label),
                            div(class = "small", style = "color: var(--text-secondary);", metric_display_label(bm$metric))),
                        info_row("Predictors", length(bm$predictors)),
                        info_row("Trees", if (!is.null(bm$reduction$n_trees_full)) {
                                sprintf("%s of %s", format(bm$n_trees, big.mark = ","), format(bm$reduction$n_trees_full, big.mark = ","))
                        } else format(bm$n_trees, big.mark = ",")),
                        if (!is.null(bm$reduction$n_trees_full))
                                div(class = "small text-muted mb-2",
                                    sprintf(paste0("Reduced forest: the first %s of the fitted forest's %s trees, unchanged. ",
                                                   "Fewer trees add Monte Carlo noise to quantiles and percentiles; the intervals include it."),
                                            format(bm$n_trees, big.mark = ","), format(bm$reduction$n_trees_full, big.mark = ",")),
                                    if (nzchar(FULL_MODELS_URL))
                                            tagList(" ", tags$a(href = FULL_MODELS_URL, target = "_blank", rel = "noopener", "Full models"),
                                                    " run in this app locally.")),
                        if (!is.null(tr))
                                info_row("Training typologies", sprintf("%s (%s data sets)", format(tr$n_rows, big.mark = ","),
                                                                        tr$n_datasets)),
                        if (!is.null(cov90) && nrow(cov90))
                                info_row("Test coverage of 90% interval", sprintf("%.0f%%", 100 * cov90$coverage[1])),
                        info_row("Calibration", if (!is.null(bm$calibration))
                                sprintf("%s held-out data sets", bm$calibration$n_datasets) else "not available"),
                        hr(),
                        div(class = "small", tags$strong(style = "color: var(--river-deep);", "Variables:"),
                            p(class = "text-muted", paste(bm$predictors, collapse = ", ")))
                )
        })
        
        # ---- Prediction ------------------------------------------------------
        collect_specified <- function(preds) {
                vals <- list()
                for (p in preds) {
                        v <- input[[paste0("param_", p)]]
                        if (is_toggled(p) && is.numeric(v) && length(v) == 1 && is.finite(v)) vals[[p]] <- v
                }
                vals
        }
        
        observeEvent(input$predict, {
                bm <- tryCatch(benchmark_full(), error = function(e) {
                        showNotification(conditionMessage(e), type = "error", duration = 8)
                        NULL
                })
                req(bm)
                specified_vals <- collect_specified(bm$predictors)
                
                tax <- intersect(SUM_CONSTRAINTS$taxonomic$vars, bm$predictors)
                if (all(tax %in% names(specified_vals)) && abs(sum(unlist(specified_vals[tax])) - 1) > 0.01) {
                        showNotification("The taxonomic rank shares you entered do not sum to 1.", type = "warning", duration = 6)
                }
                
                taus <- sort(as.numeric(input$quantiles))
                if (!length(taus)) taus <- c(0.1, 0.5, 0.9)
                obs <- input$observed_value
                obs <- if (is.numeric(obs) && length(obs) == 1 && is.finite(obs)) obs else NULL
                
                res <- withProgress(message = "Predicting", value = 0.1, tryCatch({
                        inp <- resolve_inputs(bm, specified_vals, input$impute_method, input$knn_k %||% 50)
                        incProgress(0.2, detail = "Forest prediction")
                        v <- tree_values(bm$model, inp$row[, bm$predictors, drop = FALSE])[1, ]
                        incProgress(0.2, detail = "Quantile uncertainty")
                        tg   <- tree_groups(bm, v)
                        summ <- summarise_benchmark(v, taus, obs, bm$calibration, tg$groups, tg$n_groups)
                        incProgress(0.2, detail = "Imputation effects")
                        imp <- imputation_analysis(bm, inp, v, taus, obs)
                        incProgress(0.2, detail = "Density")
                        xr  <- range(c(v, if (!is.null(imp)) range(imp$mixture$values)))
                        list(taus = taus, inputs = inp, summary = summ, imputation = imp, observed_value = obs,
                             support = tg$support, outside = outside_range(bm, specified_vals),
                             density = density_from_values(v, xr[1], xr[2]),
                             mix_density = if (!is.null(imp)) density_from_values(imp$mixture$values, xr[1], xr[2]),
                             metric = bm$metric, taxon = bm$taxon)
                }, error = function(e) {
                        showNotification(paste("Prediction error:", conditionMessage(e)), type = "error", duration = 8)
                        NULL
                }))
                prediction_result(res)
        })
        
        # ---- Prediction display ----------------------------------------------
        output$prediction_display <- renderUI({
                res <- prediction_result()
                if (is.null(res)) {
                        return(div(class = "text-center p-5",
                                   icon("chart-line", class = "fa-4x mb-3", style = "color: var(--border);"),
                                   h4(style = "color: var(--text-light); font-family: 'Instrument Serif', serif;", "No prediction yet"),
                                   p(style = "color: var(--text-light);", "Configure parameters and click 'Generate Prediction'")))
                }
                s  <- res$summary
                qt <- s$quantiles
                has_cal <- isTRUE(s$calibrated)
                
                badges <- lapply(seq_len(nrow(qt)), function(i) {
                        q <- qt$tau[i]
                        cls <- if (q < 0.4) "quantile-low" else if (q > 0.6) "quantile-high" else "quantile-mid"
                        span(class = paste("quantile-badge", cls), sprintf("Q%g: %s", q * 100, fmt(qt$value[i])))
                })
                
                ci_cell <- function(v, lo, hi) tags$td(fmt(v), tags$span(class = "uncert-ci", sprintf(" [%s, %s]", fmt(lo), fmt(hi))))
                uncert <- tagList(
                        tags$table(class = "uncert-table",
                                   tags$thead(tags$tr(tags$th("Quantile"), tags$th("Estimate"),
                                                      if (has_cal) tags$th("Calibrated"))),
                                   tags$tbody(lapply(seq_len(nrow(qt)), function(i) tags$tr(
                                           tags$td(sprintf("Q%g", 100 * qt$tau[i])),
                                           ci_cell(qt$value[i], qt$lo[i], qt$hi[i]),
                                           if (has_cal) ci_cell(qt$cal_value[i], qt$cal_lo[i], qt$cal_hi[i]))))),
                        div(class = "uncert-note",
                            sprintf("Brackets: %d%% intervals. ", level_pct(CI_LEVEL)),
                            if (isTRUE(s$grouped)) {
                                    sprintf(paste0("Estimate: each of the %s trees predicts one training value; the interval comes ",
                                                   "from resampling the trees and the training data sets those values come from ",
                                                   "(forest structure held fixed). "), format(s$n_trees, big.mark = ","))
                            } else {
                                    sprintf("Estimate: interval from resampling the forest's %s trees only (training data sets unknown). ",
                                            format(s$n_trees, big.mark = ","))
                            },
                            if (has_cal) {
                                    cal <- benchmark()$calibration
                                    sprintf(paste0("Calibrated: each quantile adjusted so that the stated share of typologies from %d ",
                                                   "held-out data sets (%s typologies) fell below it; interval from resampling those ",
                                                   "data sets and the trees."),
                                            cal$n_datasets, format(cal$n_rows, big.mark = ","))
                            } else "Calibration needs the model's _intervals.rds file, which is not available.")
                )
                
                sup <- res$support
                support_line <- if (!is.null(sup)) {
                        low <- sup$n_eff < LOW_SUPPORT
                        div(class = "small mt-2", style = paste0("color: ", if (low) "#b85c3a" else "var(--text-secondary)", ";"),
                            icon(if (low) "triangle-exclamation" else "layer-group"),
                            sprintf(" Data support: the trees' predictions come from %d training data sets (effective number %.1f; largest %.0f%%).",
                                    sup$n_datasets, sup$n_eff, 100 * sup$top_share),
                            if (low) " Few independent data sets: treat this benchmark with caution.",
                            if (sup$ambiguous > 0.2) sprintf(" %.0f%% of tree values occur in several data sets, so support is approximate.",
                                                             100 * sup$ambiguous))
                }
                
                outside_line <- if (!is.null(res$outside)) {
                        div(class = "small mt-2", style = "color: #b85c3a;", icon("triangle-exclamation"),
                            " Outside the training range (the forest does not extrapolate; predictions stay flat beyond it): ",
                            lapply(seq_len(nrow(res$outside)), function(i) {
                                    o <- res$outside[i]
                                    span(class = "imputed-var", sprintf("%s = %s [%s\u2013%s]", spec_label(o$var), fmt_in(o$value),
                                                                        fmt_in(o$min), fmt_in(o$max)))
                            }))
                }
                
                imp_chips <- NULL
                if (length(res$inputs$imputed)) {
                        method_label <- switch(res$inputs$method, conditional = "conditional averages (KNN)",
                                               marginal = "training medians", defaults = "hard-coded defaults")
                        eff <- res$imputation$effects
                        imp_chips <- div(
                                class = "imputed-summary",
                                div(class = "imputed-header", icon("sliders"),
                                    sprintf(" %d imputed via %s", length(res$inputs$imputed), method_label),
                                    if (!is.null(eff)) sprintf(" \u00b7 plausible range in brackets")),
                                div(lapply(res$inputs$imputed, function(v) {
                                        val <- res$inputs$row[[v]]
                                        rng <- if (!is.null(eff)) eff[var == v] else NULL
                                        span(class = "imputed-var", sprintf("%s = %s", spec_label(v), fmt_in(val)),
                                             if (!is.null(rng) && nrow(rng))
                                                     span(class = "imputed-range", sprintf(" [%s\u2013%s]", fmt_in(rng$lo), fmt_in(rng$hi))))
                                }))
                        )
                }
                
                div(class = "text-center py-4",
                    div(class = "prediction-value", fmt(s$median)),
                    p(style = "color: var(--text-light);", class = "mb-4", "Median Prediction"),
                    div(class = "mb-2", badges),
                    if (nrow(qt) >= 2) div(class = "small", style = "color: var(--text-light);",
                                           sprintf("%.0f%% Prediction Interval: [%s, %s]",
                                                   100 * (max(qt$tau) - min(qt$tau)), fmt(min(qt$value)), fmt(max(qt$value)))),
                    uncert,
                    support_line,
                    outside_line,
                    imp_chips)
        })
        
        output$density_subtitle <- renderUI({
                res <- prediction_result()
                txt <- "Kernel density of the forest's per-tree predictions."
                if (!is.null(res)) {
                        txt <- paste(txt, sprintf("Whiskers: %d%% interval of each quantile%s.", level_pct(CI_LEVEL),
                                                  if (isTRUE(res$summary$calibrated)) ", calibrated" else ""))
                        if (!is.null(res$mix_density)) txt <- paste(txt, "Two curves: see the explanation below the plot.")
                }
                div(class = "density-plot-subtitle", txt)
        })
        
        # Explains the solid and dashed curve whenever inputs were imputed
        output$curve_explainer <- renderUI({
                res <- prediction_result()
                imp <- res$imputation
                if (is.null(res) || is.null(imp) || is.null(res$mix_density)) return(NULL)
                lvl <- level_pct(diff(RANGE_PROBS))
                mix <- imp$mixture
                o   <- res$summary$observed
                n_imp <- length(res$inputs$imputed)
                source <- if (res$inputs$method == "conditional") "training typologies most similar to your specified inputs"
                else "randomly drawn training typologies"
                
                div(class = "curve-explainer",
                    div(class = "imputed-header", icon("circle-info"),
                        sprintf(" Why two curves? %d input%s %s imputed", n_imp, if (n_imp == 1) "" else "s",
                                if (n_imp == 1) "was" else "were")),
                    div(class = "curve-row",
                        span(class = "curve-key"),
                        div(tags$b(paste0(CURVE_LABELS[["fixed"]], ":")), " each imputed input takes one value, as if it were known. ",
                            sprintf("Median %s, %d%% interval [%s, %s]. ", fmt(imp$base_median), lvl,
                                    fmt(imp$base_range[1]), fmt(imp$base_range[2])),
                            "The quantile lines, whiskers and percentile shown here refer to this curve.")),
                    div(class = "curve-row",
                        span(class = "curve-key dashed"),
                        div(tags$b(paste0(CURVE_LABELS[["varied"]], ":")),
                            sprintf(" the imputed inputs take their values from %d %s, so the curve shows the benchmark you can expect when they are not known. ",
                                    imp$n_draws, source),
                            paste0(sprintf("Median %s, %d%% interval [%s, %s]", fmt(mix$median), lvl, fmt(mix$range[1]), fmt(mix$range[2])),
                                   if (!is.null(o) && !is.null(mix$pct)) sprintf("; your observed value is at its %s percentile.", pct_label(mix$pct)) else "."))),
                    div(class = "curve-reading",
                        tags$b("Reading them: "),
                        "a wider dashed curve means the unknown inputs add uncertainty; a shifted one means the single best guesses ",
                        "do not represent them well; overlapping curves mean the imputed inputs hardly matter. ",
                        "Specifying more inputs brings the two curves together. ",
                        "The card \u201cImputed Inputs and Their Effect\u201d shows which inputs matter most."))
        })
        
        prediction_plot <- function(res) {
                qt <- res$summary$quantiles
                qtab <- if (isTRUE(res$summary$calibrated)) {
                        data.frame(tau = qt$tau, value = qt$value, lo = qt$cal_lo, hi = qt$cal_hi)
                } else {
                        data.frame(tau = qt$tau, value = qt$value, lo = qt$lo, hi = qt$hi)
                }
                build_density_plot(res$density, qtab, res$observed_value, res$mix_density,
                                   metric_display_label(res$metric))
        }
        
        output$density_plot <- renderPlot({
                res <- prediction_result()
                if (is.null(res) || is.null(res$density)) {
                        return(placeholder_plot("Generate a prediction to see the distribution"))
                }
                prediction_plot(res)
        }, bg = "transparent", res = 96)
        
        # Downloaded figure: the plot plus title, key numbers and an explanation
        output$download_density <- downloadHandler(
                filename = function() {
                        res <- prediction_result()
                        stem <- if (is.null(res)) "prediction" else paste0(res$taxon, "__", res$metric)
                        paste0(stem, "_distribution_", Sys.Date(), ".", input$figure_format %||% "png")
                },
                content = function(file) {
                        res <- prediction_result()
                        validate(need(!is.null(res) && !is.null(res$density), "Generate a prediction first."))
                        s   <- res$summary
                        lvl <- level_pct(diff(RANGE_PROBS))
                        sub <- sprintf("Median %s, %d%% interval [%s, %s]", fmt(s$median), lvl, fmt(s$range[1]), fmt(s$range[2]))
                        if (!is.null(s$observed)) {
                                o <- s$observed
                                sub <- paste0(sub, sprintf("; observed value %s at the %s percentile", fmt(o$value),
                                                           pct_label(o$cal_pct %||% o$pct)),
                                              if (!is.null(o$cal_pct)) " (calibrated)" else "")
                        }
                        n_imp <- length(res$inputs$imputed)
                        cap <- paste0(
                                sprintf("Kernel density of the forest's per-tree predictions. Vertical lines: quantiles; whiskers: %d%% interval of each quantile%s. ",
                                        level_pct(CI_LEVEL), if (isTRUE(s$calibrated)) " (calibrated to held-out data sets)" else ""),
                                if (n_imp) sprintf("%d of %d inputs imputed (%s)%s. ", n_imp, n_imp + length(res$inputs$specified),
                                                   switch(res$inputs$method, conditional = "KNN", marginal = "training medians",
                                                          defaults = "defaults"),
                                                   if (!is.null(res$mix_density)) "; dashed: imputed inputs varied over plausible values" else "")
                                else "All inputs specified. ",
                                "PULSE benchmark model, ", format(Sys.Date()), ".")
                        p <- prediction_plot(res) +
                                labs(title = sprintf("%s - %s", taxon_info(res$taxon)$label, metric_display_label(res$metric)),
                                     subtitle = sub, caption = wrap_text(cap, 120)) +
                                theme(plot.title = element_text(face = "bold", size = 13),
                                      plot.subtitle = element_text(size = 10, color = PAL$text),
                                      plot.caption = element_text(size = 8, color = PAL$light, hjust = 0))
                        save_plot(p, file, input$figure_format %||% "png", width = 9, height = 5.8)
                }
        )
        
        output$observed_percentile_display <- renderUI({
                res <- prediction_result()
                o <- res$summary$observed
                if (is.null(o)) return(NULL)
                head_pct  <- o$cal_pct %||% o$pct
                pct       <- 100 * head_pct
                pct_color <- if (pct < 10 || pct > 90) "#b85c3a" else if (pct < 25 || pct > 75) "#c4a882" else "#1e7a5a"
                div(class = "text-center mt-2 mb-3",
                    div(class = "percentile-badge", style = paste0("background: ", pct_color, "15; color: ", pct_color, ";"),
                        icon("crosshairs"),
                        sprintf(" Observed value (%s) is at the %s percentile [%s]", fmt(o$value), pct_label(o$pct), pct_range(o$pct_ci))),
                    if (!is.null(o$cal_pct))
                            div(class = "small mt-1", style = "color: var(--text-secondary);",
                                sprintf("Calibrated to held-out data sets: %s percentile [%s]", pct_label(o$cal_pct), pct_range(o$cal_ci))),
                    div(class = "small mt-1", style = paste0("color: ", pct_color, ";"),
                        if (pct < 10) "Below the 10th percentile \u2014 substantially lower than expected"
                        else if (pct < 25) "Below the 25th percentile \u2014 lower than typical"
                        else if (pct > 90) "Above the 90th percentile \u2014 substantially higher than expected"
                        else if (pct > 75) "Above the 75th percentile \u2014 higher than typical"
                        else "Within the interquartile range \u2014 consistent with expectations"))
        })
        
        # ---- Imputation panel ------------------------------------------------
        output$imputation_panel <- renderUI({
                res <- prediction_result()
                if (is.null(res)) {
                        return(p(class = "text-muted p-2", "Generate a prediction to see which inputs were imputed and how much they affect the benchmark."))
                }
                inp <- res$inputs
                if (!length(inp$imputed)) {
                        return(p(class = "p-2", style = "color: var(--accent);", icon("circle-check"),
                                 " All model inputs were specified, so imputation adds no uncertainty."))
                }
                imp <- res$imputation
                if (is.null(imp)) {
                        return(p(class = "p-2", style = "color: var(--text-secondary);", icon("triangle-exclamation"),
                                 sprintf(" %d inputs were set to hard-coded defaults. Ranges and effects need the training data, which this model file does not carry.",
                                         length(inp$imputed))))
                }
                lvl    <- level_pct(diff(RANGE_PROBS))
                source <- if (inp$method == "conditional") sprintf("the %d training typologies nearest to your specified inputs", imp$pool_size)
                else sprintf("%s training typologies", format(imp$pool_size, big.mark = ","))
                mix    <- imp$mixture
                n_show <- min(N_TORNADO, nrow(imp$effects))
                
                tagList(
                        div(class = "imputation-text",
                            span(class = paste("impact-badge", paste0("impact-", imp$impact)), paste("Impact of imputation:", imp$impact)),
                            p(sprintf("%d of %d inputs were imputed. Their plausible ranges are the 10th\u201390th percentiles among %s.",
                                      length(inp$imputed), length(inp$imputed) + length(inp$specified), source)),
                            p("Across ", imp$n_draws, " plausible combinations of the imputed inputs, the median benchmark lies between ",
                              tags$b(fmt(imp$band[1])), " and ", tags$b(fmt(imp$band[2])),
                              sprintf(" (%d%% of combinations). This spread is ", lvl),
                              tags$b(if (is.na(imp$ratio)) "NA" else sprintf("%.0f%%", 100 * imp$ratio)),
                              sprintf(" of the width of the benchmark's own %d%% interval [%s, %s].", lvl,
                                      fmt(imp$base_range[1]), fmt(imp$base_range[2]))),
                            p(paste0(sprintf("Benchmark accounting for imputation: median %s, %d%% interval [%s, %s]",
                                             fmt(mix$median), lvl, fmt(mix$range[1]), fmt(mix$range[2])),
                                     if (!is.null(mix$pct)) sprintf("; your observed value is at its %s percentile.", pct_label(mix$pct)) else "."))),
                        div(class = "density-plot-header mt-3", "Effect of each imputed input"),
                        div(class = "density-plot-subtitle",
                            "Each bar shows how the median benchmark changes when one input moves from the lower (\u25cb) to the upper (\u25cf) ",
                            "end of its plausible range, the other inputs held at their imputed values. ",
                            "Shaded band (dashed edges): spread of the median over plausible combinations of all imputed inputs together."),
                        plotOutput("tornado_plot", height = paste0(tornado_height(n_show), "px"), fill = FALSE),
                        if (nrow(imp$effects) > n_show)
                                tags$small(class = "text-muted", sprintf("%d further imputed inputs with smaller effects are not shown.",
                                                                         nrow(imp$effects) - n_show))
                )
        })
        
        output$tornado_plot <- renderPlot({
                res <- prediction_result()
                req(res$imputation)
                build_tornado_plot(res$imputation, metric_display_label(res$metric))
        }, height = function() {
                eff <- prediction_result()$imputation$effects
                tornado_height(if (is.null(eff)) 1 else min(N_TORNADO, nrow(eff)))
        }, bg = "transparent", res = 96)
        
        # ---- Batch -----------------------------------------------------------
        batch_data <- reactive({
                req(input$batch_file)
                tryCatch(fread(input$batch_file$datapath), error = function(e) NULL)
        })
        
        output$batch_column_status <- renderUI({
                df <- batch_data()
                if (is.null(df)) return(div(class = "small mt-2", style = "color: var(--danger);", "The file could not be read as CSV."))
                req(input$batch_taxon_group)
                res  <- resolve_batch_models(as.data.frame(df), input$batch_taxon_group, input$batch_metric)
                pred <- names(PREDICTOR_SPEC)
                found <- intersect(pred, names(df)); missing <- setdiff(pred, names(df))
                ok   <- res$status == "ok"
                tab  <- if (any(ok)) table(sprintf("%s - %s", vapply(res$taxon[ok], function(t) taxon_info(t)$label, ""),
                                                   metric_display_label(res$metric[ok]))) else NULL
                div(class = "imputed-summary mt-3",
                    div(class = "imputed-header", icon("table-columns"),
                        sprintf(" %d rows | %d / %d predictors found", nrow(df), length(found), length(pred))),
                    if (length(missing))
                            div(div(class = "small mt-1 mb-1", style = "color: var(--water-blue);",
                                    sprintf("%d missing \u2014 will be imputed:", length(missing))),
                                div(lapply(missing, function(v) span(class = "imputed-var", spec_label(v))))),
                    div(class = "small mt-2", style = "color: var(--text-secondary);",
                        tags$b(if (is.null(res$metric_col)) "Model (from the selection above): " else "Models (from the metric column): "),
                        if (length(tab)) paste(sprintf("%s: %d rows", names(tab), as.integer(tab)), collapse = "; ") else "none"),
                    if (any(!ok))
                            div(class = "small mt-1", style = "color: var(--danger);", icon("triangle-exclamation"),
                                sprintf(" %d rows without a model (%s); they are kept in the results without predictions.",
                                        sum(!ok), paste(unique(res$status[!ok]), collapse = ", "))),
                    if (!is.null(find_col(df, "observed_value")))
                            div(class = "small mt-1", style = "color: var(--accent);", icon("crosshairs"),
                                " observed_value found: percentiles will be added"))
        })
        
        observeEvent(input$batch_predict, {
                df <- batch_data()
                req(df)
                taus <- sort(as.numeric(input$batch_quantiles))
                if (!length(taus)) taus <- c(0.1, 0.5, 0.9)
                withProgress(message = "Batch prediction", value = 0, tryCatch({
                        b <- run_batch(df, input$batch_taxon_group, input$batch_metric, taus,
                                       input$batch_impute_method, input$batch_knn_k %||% 50,
                                       progress = function(amount, detail = NULL) incProgress(0.9 * amount, detail = detail))
                        batch_results(b)
                        n_ok <- sum(b$table$status == "ok")
                        showNotification(sprintf("Predictions complete: %d of %d rows, %d model%s", n_ok, nrow(b$table),
                                                 length(b$groups), if (length(b$groups) == 1) "" else "s"),
                                         type = if (n_ok < nrow(b$table)) "warning" else "message", duration = 6)
                }, error = function(e) {
                        showNotification(paste("Error:", conditionMessage(e)), type = "error", duration = 8)
                }))
        })
        
        output$batch_results <- renderDT({
                req(batch_results())
                tb <- batch_results()$table
                first <- c(intersect(c("typology", "id", "name", "label", "site"), names(tb))[1],
                           "model_taxon", "model_metric", "status", grep("^Q[0-9]", names(tb), value = TRUE),
                           intersect(c("observed_value", "observed_percentile", "observed_percentile_calibrated",
                                       "support_datasets_eff", "imputed_inputs"), names(tb)))
                first <- first[!is.na(first)]
                tb <- tb[, c(first, setdiff(names(tb), first)), with = FALSE]
                num <- names(tb)[vapply(tb, is.double, TRUE)]
                datatable(tb, options = list(pageLength = 10, scrollX = TRUE), class = "cell-border stripe",
                          rownames = FALSE, fillContainer = FALSE) |> formatSignif(num, digits = 4)
        })
        
        output$batch_figure_panel <- renderUI({
                b <- batch_results()
                if (is.null(b) || !length(b$groups)) {
                        return(p(class = "text-muted p-2", "Run a batch prediction to see an overview figure per metric."))
                }
                choices <- setNames(names(b$groups), vapply(b$groups, group_title, ""))
                tagList(
                        div(class = "d-flex flex-wrap align-items-end gap-3",
                            selectInput("batch_figure_group", "Model", choices = choices, width = "320px"),
                            radioButtons("batch_figure_format", NULL, choices = FIGURE_FORMATS, selected = "png", inline = TRUE),
                            downloadButton("download_batch_figure", "Download this figure", class = "btn-outline-secondary btn-sm mb-3")),
                        uiOutput("batch_overview_container")
                )
        })
        
        batch_group <- reactive({
                b <- batch_results(); req(b, input$batch_figure_group)
                b$groups[[input$batch_figure_group]]
        })
        
        output$batch_overview_container <- renderUI({
                g <- batch_group(); req(g)
                plotOutput("batch_overview", height = paste0(round(96 * batch_overview_height(g)), "px"), fill = FALSE)
        })
        
        output$batch_overview <- renderPlot({
                g <- batch_group(); req(g)
                build_batch_overview(g, title = FALSE)
        }, height = function() round(96 * batch_overview_height(batch_group())), bg = "transparent", res = 96)
        
        output$download_batch_figure <- downloadHandler(
                filename = function() paste0(input$batch_figure_group, "_overview.", input$batch_figure_format %||% "png"),
                content  = function(file) {
                        g <- batch_group()
                        save_plot(build_batch_overview(g), file, input$batch_figure_format %||% "png",
                                  width = 8.5, height = batch_overview_height(g))
                }
        )
        
        output$download_batch_archive <- downloadHandler(
                filename = function() paste0("pulse_batch_", Sys.Date(), ".", ARCHIVE_EXT),
                content  = function(file) {
                        b <- batch_results()
                        validate(need(!is.null(b), "Run a batch prediction first."))
                        dir <- file.path(tempfile("pulse_batch_"), "pulse_batch")
                        withProgress(message = "Writing files and figures", value = 0.3, write_batch_outputs(b, dir))
                        make_archive(dir, file)
                }
        )
        
        output$download_results <- downloadHandler(
                filename = function() paste0("pulse_batch_all_results_", Sys.Date(), ".csv"),
                content  = function(file) { req(batch_results()); fwrite(batch_results()$table, file) }
        )
        
        # ---- Parameter reference ---------------------------------------------
        output$param_reference <- renderDT({
                ref_data <- data.frame(
                        Parameter   = names(PREDICTOR_SPEC),
                        Label       = spec_label(names(PREDICTOR_SPEC)),
                        Group       = vapply(PREDICTOR_SPEC, `[[`, "", "group"),
                        Description = vapply(PREDICTOR_SPEC, function(x) x$help %||% "", ""),
                        Default     = vapply(PREDICTOR_SPEC, function(x) as.character(x$value), "")
                )
                datatable(ref_data, options = list(pageLength = 25, dom = "t", ordering = FALSE),
                          rownames = FALSE, class = "compact")
        })
}


# Run App -----------------------------------------------------------------
shinyApp(ui = ui, server = server)