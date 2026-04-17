### Quantile Random Forest Prediction App ###


# 

# Startup Diagnostics -----------------------------------------------------

cat("\n\n")
cat("QRF Prediction App - Starting...\n")
cat("\n\n")

# Check and load packages with informative errors
required_packages <- c(
        "shiny", "bslib", "workflows", "recipes", "parsnip", "bundle", 
        "ranger", "data.table", "DT", "ggplot2"
)

missing_packages <- c()

for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
                missing_packages <- c(missing_packages, pkg)
                cat("MISSING:", pkg, "\n")
        } else {
                cat("OK:", pkg, "\n")
        }
}

if (length(missing_packages) > 0) {
        stop(
                "Missing required packages: ", 
                paste(missing_packages, collapse = ", "),
                "\nInstall with: install.packages(c('", 
                paste(missing_packages, collapse = "', '"), "'))"
        )
}

cat("\nAll packages available. Loading...\n\n")

library(shiny)
library(bslib)
library(workflows)
library(recipes)
library(parsnip)
library(bundle)
library(ranger)
library(data.table)
library(DT)
library(ggplot2)

# Increase time before time-out error when downloading files. 
options(timeout = max(600, getOption("timeout")))

# Configuration -----------------------------------------------------------

# Base URL for the Zenodo record. All model files live directly under /files/
ZENODO_BASE_URL <- "https://zenodo.org/records/19554225/files"
# Local cache directory — persists across sessions on the same host
#CACHE_DIR <- tools::R_user_dir("qrf_predictor", which = "cache")
# SERVER VERSION
CACHE_DIR <- Sys.getenv("CACHE_DIR", "/srv/shiny-server/cache")
if (!dir.exists(CACHE_DIR)) {
        dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
}
cat("Cache directory:", CACHE_DIR, "\n")


# MODEL MANIFEST ----------------------------------------------------------
# Manifest: which metrics exist for which taxonomic group, and whether a
# shared training_data file is available. You maintain this by hand to match
# what you've uploaded to Zenodo. The app no longer auto-discovers files
# because there's no directory to scan.
#
# File naming convention on Zenodo (flat layout, no subfolders):
#   <taxon>__<metric>.rds            e.g. macroinvertebrates__ANOSIM_R_max.rds
#   <taxon>__training_data.rds       e.g. macroinvertebrates__training_data.rds
#
# If your Zenodo deposit uses a different naming scheme, adjust
# build_remote_filename() below to match.
MODEL_MANIFEST <- list(
        diatoms = c(
                "ANOSIM_R_max",
                "ANOSIM_R_mean",
                "ANOSIM_R_min",
                "AucZeta_max",
                "AucZeta_mean",
                "AucZeta_min",
                "classification_strength",
                "fuzzy_mantel",
                "isamic",
                "isa_avg_p",
                "isa_number",
                "PERMANOVA_Fuzzy_R2",
                "PERMANOVA_R2"
        ),
        fishes = c(
                "ANOSIM_R_max",
                "ANOSIM_R_mean",
                "ANOSIM_R_min",
                "AucZeta_max",
                "AucZeta_mean",
                "AucZeta_min",
                "classification_strength",
                "fuzzy_mantel",
                "isamic",
                "isa_avg_p",
                "isa_number",
                "PERMANOVA_Fuzzy_R2",
                "PERMANOVA_R2"
        ),
        macroinvertebrates = c(
                "ANOSIM_R_max",
                "ANOSIM_R_mean",
                "ANOSIM_R_min",
                "AucZeta_max",
                "AucZeta_mean",
                "AucZeta_min",
                "classification_strength",
                "fuzzy_mantel",
                "isamic",
                "isa_avg_p",
                "isa_number",
                "PERMANOVA_Fuzzy_R2",
                "PERMANOVA_R2"
        ),
        macrophytes = c(
                "ANOSIM_R_max",
                "ANOSIM_R_mean",
                "ANOSIM_R_min",
                "AucZeta_max",
                "AucZeta_mean",
                "AucZeta_min",
                "classification_strength",
                "fuzzy_mantel",
                "isamic",
                "isa_avg_p",
                "isa_number",
                "PERMANOVA_Fuzzy_R2",
                "PERMANOVA_R2"
        )
)


# Remove cache functions --------------------------------------------------

#' Remove all cached files. Returns a summary of what was removed.
clear_cache <- function() {
        if (!dir.exists(CACHE_DIR)) {
                return(list(n_files = 0, total_mb = 0))
        }
        files <- list.files(CACHE_DIR, full.names = TRUE)
        if (length(files) == 0) {
                return(list(n_files = 0, total_mb = 0))
        }
        total_bytes <- sum(file.size(files), na.rm = TRUE)
        file.remove(files)
        list(n_files = length(files), total_mb = total_bytes / 1e6)
}

#' Summarise current cache contents without clearing.
cache_summary <- function() {
        if (!dir.exists(CACHE_DIR)) {
                return(list(n_files = 0, total_mb = 0))
        }
        files <- list.files(CACHE_DIR, full.names = TRUE)
        list(
                n_files = length(files),
                total_mb = if (length(files) > 0) sum(file.size(files), na.rm = TRUE) / 1e6 else 0
        )
}



# Remote file fetching with local cache


#' Build the remote filename for a given taxon + metric
#' Adjust this to match your actual Zenodo file naming convention.
build_remote_filename <- function(taxon_group, metric = NULL, 
                                  type = c("model", "training_data")) {
        type <- match.arg(type)
        if (type == "model") {
                sprintf("%s_%s.rds", taxon_group, metric)
        } else {
                sprintf("%s_training_data.rds", taxon_group)
        }
}

#' Download a file from Zenodo into the local cache (if not already present)
#' and return the local path. Returns NULL on failure.
fetch_remote_file <- function(filename, force = FALSE) {
        local_path <- file.path(CACHE_DIR, filename)
        
        if (file.exists(local_path) && !force) {
                cat("  Cache hit:", filename, "\n")
                return(local_path)
        }
        
        url <- paste0(ZENODO_BASE_URL, "/", filename)
        cat("  Downloading:", url, "\n")
        
        tryCatch({
                # mode = "wb" is critical for binary .rds files on Windows
                download.file(url, destfile = local_path, mode = "wb", quiet = TRUE)
                if (file.exists(local_path) && file.size(local_path) > 0) {
                        cat("  Downloaded:", filename, 
                            sprintf("(%.1f MB)", file.size(local_path) / 1e6), "\n")
                        return(local_path)
                } else {
                        return(NULL)
                }
        }, error = function(e) {
                cat("  Download failed:", e$message, "\n")
                if (file.exists(local_path)) file.remove(local_path)  # don't leave partials
                return(NULL)
        })
}

#' Get available metrics for a taxon from the manifest (replaces filesystem scan)
get_available_metrics <- function(taxon_group = NULL) {
        if (is.null(taxon_group) || taxon_group == "" || 
            is.null(MODEL_MANIFEST[[taxon_group]])) {
                return(character(0))
        }
        MODEL_MANIFEST[[taxon_group]]
}
cat("Model source: Zenodo (", ZENODO_BASE_URL, ")\n")
cat("Cache location:", CACHE_DIR, "\n\n")

# Taxonomic groups with display names and directory suffixes
TAXONOMIC_GROUPS <- list(
        diatoms = list(
                label = "Diatoms",
                dir = "diatoms",
                icon = "leaf"
        ),
        fishes = list(
                label = "Fishes",
                dir = "fishes",
                icon = "fish"
        ),
        macroinvertebrates = list(
                label = "Macroinvertebrates",
                dir = "macroinvertebrates",
                icon = "bug"
        ),
        macrophytes = list(
                label = "Macrophytes",
                dir = "macrophytes",
                icon = "seedling"
        )
)

# Define groups of variables that must sum to 1
COMPOSITIONAL_GROUPS <- list(
        variation = c("bio", "env", "space", "stochastic"),
        taxonomic = c("species_rank", "genus_rank", "family_rank", "higher_rank")
)

# Ordered variable chains: each imputed value is clamped to its neighbors
ORDERED_CONSTRAINTS <- list(
        distance = list(
                chain = c("min_distance", "median_distance", "mean_distance", "max_distance"),
                interior = c("median_distance", "mean_distance"),
                endpoints = c("min_distance", "max_distance")
        ),
        latitude = list(
                chain = c("min_latitude", "median_latitude", "max_latitude"),
                interior = c("median_latitude"),
                endpoints = c("min_latitude", "max_latitude")
        ),
        longitude = list(
                chain = c("min_longitude", "median_longitude", "max_longitude"),
                interior = c("median_longitude"),
                endpoints = c("min_longitude", "max_longitude")
        )
)

# Conditional Averaging via KNN -------------------------------------------

#' Compute conditional averages for unspecified variables using KNN
#'
#' Given a training dataset, specified variable values, and which variables
#' are unspecified, find the k nearest neighbors in the training data based
#' on the specified variables and return the average of the unspecified 
#' variables from those neighbors.
#'
#' @param training_data data.frame of training data (predictors only)
#' @param specified_vals named list of specified variable values
#' @param unspecified_vars character vector of variable names to impute
#' @param k number of nearest neighbors (default: 50, capped at nrow/5)
#' @return named numeric vector of conditional averages for unspecified vars
compute_conditional_averages <- function(training_data, specified_vals, 
                                         unspecified_vars, k = 50) {
        
        specified_vars <- names(specified_vals)
        
        # Ensure all specified vars exist in training data
        specified_vars <- intersect(specified_vars, names(training_data))
        unspecified_vars <- intersect(unspecified_vars, names(training_data))
        
        if (length(specified_vars) == 0 || nrow(training_data) == 0) {
                # Fall back to marginal medians
                return(sapply(unspecified_vars, function(v) {
                        median(training_data[[v]], na.rm = TRUE)
                }))
        }
        
        if (length(unspecified_vars) == 0) {
                return(numeric(0))
        }
        
        # Cap k at a reasonable fraction of training data
        k <- min(k, max(10, floor(nrow(training_data) / 5)))
        k <- min(k, nrow(training_data))
        
        # Extract the specified columns from training data
        train_spec <- as.matrix(training_data[, specified_vars, drop = FALSE])
        
        # Scale using training data ranges for fair distance computation
        col_sds <- apply(train_spec, 2, sd, na.rm = TRUE)
        col_means <- colMeans(train_spec, na.rm = TRUE)
        
        # Avoid division by zero
        col_sds[col_sds == 0] <- 1
        
        # Scale training data
        train_scaled <- scale(train_spec, center = col_means, scale = col_sds)
        
        # Scale the query point
        query_vals <- unlist(specified_vals[specified_vars])
        query_scaled <- (query_vals - col_means) / col_sds
        
        # Compute Euclidean distances
        diffs <- sweep(train_scaled, 2, query_scaled, "-")
        dists <- sqrt(rowSums(diffs^2, na.rm = TRUE))
        
        # Find k nearest neighbors
        knn_idx <- order(dists)[1:k]
        
        # Compute averages of unspecified variables from neighbors
        neighbor_data <- training_data[knn_idx, unspecified_vars, drop = FALSE]
        
        cond_avgs <- colMeans(neighbor_data, na.rm = TRUE)
        names(cond_avgs) <- unspecified_vars
        
        return(cond_avgs)
}




# estimate_density_from_qrf() ----------------------------

#' Estimate a density curve from a dense quantile prediction grid
#'
#' Predicts on a fine quantile grid (e.g., 0.01 to 0.99) using the ranger
#' QRF model, then numerically differentiates the inverse-CDF to obtain
#' an approximate density function.
#'
#' @param ranger_model Fitted ranger model (extracted from workflow)
#' @param prepped_data Baked/prepped single-row data.frame for prediction
#' @param from Lower quantile bound (default 0.005)
#' @param to Upper quantile bound (default 0.995)
#' @param n_quantiles Number of quantile evaluation points (default 199)
#' @param smooth Whether to apply loess smoothing (default TRUE)
#' @param smooth_span Loess span parameter (default 0.2)
#' @return data.frame with columns: x (metric value), density, quantile
# estimate_density_from_qrf <- function(ranger_model, prepped_data,
#                                       from = 0.005, to = 0.995,
#                                       n_quantiles = 199,
#                                       smooth = TRUE,
#                                       smooth_span = 0.35) {
#         
#         # Dense quantile grid
#         fine_quantiles <- seq(from, to, length.out = n_quantiles)
#         
#         # Predict on the fine grid
#         qpred_dense <- predict(
#                 ranger_model,
#                 data = prepped_data,
#                 type = "quantiles",
#                 quantiles = fine_quantiles
#         )
#         
#         qvals <- as.numeric(qpred_dense$predictions[1, ])
#         
#         # Numerical differentiation: density ≈ Δp / Δx
#         dx <- diff(qvals)
#         dp <- diff(fine_quantiles)
#         
#         # Avoid division by zero for flat segments
#         density_est <- ifelse(abs(dx) < .Machine$double.eps * 100, 0, dp / dx)
#         
#         # Midpoints for plotting
#         midpoints <- (qvals[-1] + qvals[-length(qvals)]) / 2
#         mid_quantiles <- (fine_quantiles[-1] + fine_quantiles[-length(fine_quantiles)]) / 2
#         
#         dens_df <- data.frame(
#                 x = midpoints,
#                 density = density_est,
#                 quantile = mid_quantiles
#         )
#         
#         # Remove non-finite or negative density values
#         dens_df <- dens_df[is.finite(dens_df$density) & dens_df$density >= 0, ]
#         
#         if (nrow(dens_df) < 5) {
#                 return(dens_df)
#         }
#         
#         # Smoothing pipeline: monotone spline on quantile function, then loess on density
#         if (smooth && nrow(dens_df) > 10) {
#                 tryCatch({
#                         # Step 1: Smooth the quantile function with a monotone spline
#                         #         before differentiating — this removes the jaggedness
#                         #         at source rather than patching it after the fact.
#                         qfun <- splinefun(fine_quantiles, qvals, method = "monoH.FC")
#                         
#                         # Evaluate on a denser, evenly-spaced grid
#                         eval_p <- seq(from + 0.005, to - 0.005, length.out = 300)
#                         eval_x <- qfun(eval_p)
#                         
#                         # Numerical derivative of the spline
#                         dx2 <- diff(eval_x)
#                         dp2 <- diff(eval_p)
#                         dens2 <- ifelse(abs(dx2) < .Machine$double.eps * 100, 0, dp2 / dx2)
#                         mid_x2 <- (eval_x[-1] + eval_x[-length(eval_x)]) / 2
#                         mid_p2 <- (eval_p[-1] + eval_p[-length(eval_p)]) / 2
#                         
#                         dens_df <- data.frame(
#                                 x = mid_x2, density = dens2, quantile = mid_p2
#                         )
#                         dens_df <- dens_df[is.finite(dens_df$density) & dens_df$density >= 0, ]
#                         
#                         # Step 2: Light loess pass for final visual polish
#                         lo <- loess(density ~ x, data = dens_df,
#                                     span = smooth_span,
#                                     degree = 2)
#                         dens_df$density_raw <- dens_df$density
#                         dens_df$density <- pmax(0, predict(lo, dens_df$x))
#                         
#                 }, error = function(e) {
#                         # If spline approach fails, fall back to loess only
#                         tryCatch({
#                                 lo <- loess(density ~ x, data = dens_df,
#                                             span = smooth_span,
#                                             degree = 2)
#                                 dens_df$density_raw <- dens_df$density
#                                 dens_df$density <- pmax(0, predict(lo, dens_df$x))
#                         }, error = function(e2) {
#                                 dens_df$density_raw <- dens_df$density
#                         })
#                 })
#         }
#         
#         return(dens_df)
# }
estimate_density_from_qrf <- function(ranger_model, prepped_data,
                                      n_samples = 500,
                                      bw_method = "nrd0") {
        
        # Predict on a dense quantile grid — this gives us an empirical
        # sample of the predictive distribution
        fine_quantiles <- seq(0.005, 0.995, length.out = n_samples)
        qpred <- predict(ranger_model, data = prepped_data,
                         type = "quantiles", quantiles = fine_quantiles)
        qvals <- as.numeric(qpred$predictions[1, ])
        
        # KDE on the sampled values. Many will be duplicates (step function),
        # but that's fine — density() handles ties correctly by design.
        # Using adjust slightly > 1 helps when the support has many ties.
        dens <- density(qvals, bw = bw_method, adjust = 1.2, n = 512,
                        from = min(qvals), to = max(qvals))
        
        data.frame(
                x = dens$x,
                density = dens$y,
                quantile = NA_real_  # no longer meaningful per point
        )
}

# build_density_plot ------------------------------------------------------

#' Build a styled density plot with optional observed value marker
#'
#' @param dens_df data.frame from estimate_density_from_qrf()
#' @param highlight_quantiles Numeric vector of quantiles to mark (e.g., c(0.1, 0.5, 0.9))
#' @param highlight_values Numeric vector of corresponding predicted values
#' @param observed_value Optional numeric observed value to mark on the plot
#' @param metric_label Label for the x-axis
#' @return ggplot object
#' 
build_density_plot <- function(dens_df,
                               highlight_quantiles = NULL,
                               highlight_values = NULL,
                               observed_value = NULL,
                               metric_label = "Metric Value") {
        
        if (is.null(dens_df) || nrow(dens_df) < 3) {
                # Return a placeholder plot
                return(
                        ggplot() +
                                annotate("text", x = 0.5, y = 0.5,
                                         label = "Insufficient data for density estimate",
                                         size = 4, color = "#7a7a7a") +
                                theme_void()
                )
        }
        
        # Colour palette matching the app theme
        col_fill   <- "#2e6b8a"
        col_line   <- "#1e7a5a"
        col_obs    <- "#b85c3a"
        col_quant  <- "#4a4a4a"
        col_bg     <- "#f5f1ea"
        col_border <- "#e0dcd4"
        col_text   <- "#4a4a4a"
        col_light  <- "#7a7a7a"
        
        p <- ggplot(dens_df, aes(x = x, y = density)) +
                geom_area(alpha = 0.18, fill = col_fill) +
                #geom_line(aes(y = density_raw), color = "red", linewidth = 0.5, alpha = 0.4) +
                geom_line(color = col_line, linewidth = 0.9)
        
        # Shade inter-quantile region if we have highlight quantiles
        if (!is.null(highlight_quantiles) && length(highlight_values) >= 2) {
                q_lo <- min(highlight_values)
                q_hi <- max(highlight_values)
                shade_df <- dens_df[dens_df$x >= q_lo & dens_df$x <= q_hi, ]
                if (nrow(shade_df) > 1) {
                        p <- p +
                                geom_area(data = shade_df,
                                          aes(x = x, y = density),
                                          alpha = 0.12, fill = col_line)
                }
        }
        
        # Add quantile markers as full-height lines
        if (!is.null(highlight_quantiles) && !is.null(highlight_values)) {
                max_dens <- max(dens_df$density, na.rm = TRUE)
                
                for (i in seq_along(highlight_quantiles)) {
                        q <- highlight_quantiles[i]
                        v <- highlight_values[i]
                        
                        is_median <- abs(q - 0.5) < 0.01
                        
                        # Lines span full plot height so they're visible even in the tails
                        p <- p +
                                geom_vline(
                                        xintercept = v,
                                        linetype = if (is_median) "solid" else "dashed",
                                        color = if (is_median) col_line else col_quant,
                                        linewidth = if (is_median) 1.0 else 0.65,
                                        alpha = if (is_median) 0.85 else 0.6
                                ) +
                                annotate(
                                        "label",
                                        x = v,
                                        y = max_dens * if (is_median) 0.97 else 0.88,
                                        label = sprintf("Q%g", q * 100),
                                        size = if (is_median) 3.4 else 2.9,
                                        color = if (is_median) col_line else col_quant,
                                        fill = col_bg,
                                        fontface = if (is_median) "bold" else "plain",
                                        label.padding = unit(0.2, "lines"),
                                        label.r = unit(0.2, "lines"),
                                        label.size = 0.2
                                )
                }
        }
        
        # Add observed value marker
        if (!is.null(observed_value) && is.finite(observed_value)) {
                # Interpolate density at observed value
                y_at_obs <- tryCatch({
                        approx(dens_df$x, dens_df$density, xout = observed_value)$y
                }, error = function(e) 0)
                if (is.null(y_at_obs) || is.na(y_at_obs)) y_at_obs <- max(dens_df$density) * 0.7
                
                p <- p +
                        geom_vline(
                                xintercept = observed_value,
                                linetype = "solid",
                                color = col_obs,
                                linewidth = 1.2
                        ) +
                        annotate(
                                "label",
                                x = observed_value,
                                y = max(dens_df$density) * 1.08,
                                label = sprintf("Observed: %.3f", observed_value),
                                fill = col_bg,
                                color = col_obs,
                                fontface = "bold",
                                size = 3.2,
                                label.padding = unit(0.35, "lines"),
                                label.r = unit(0.3, "lines"),
                                label.size = 0.3
                        )
        }
        
        # Theme
        p <- p +
                scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
                coord_cartesian(clip = "off") +
                labs(
                        x = metric_label,
                        y = "Density"
                ) +
                theme_minimal(base_family = "sans", base_size = 12) +
                theme(
                        plot.background  = element_rect(fill = "transparent", color = NA),
                        panel.background = element_rect(fill = "transparent", color = NA),
                        panel.grid.major = element_line(color = col_border, linewidth = 0.3),
                        panel.grid.minor = element_blank(),
                        axis.text        = element_text(color = col_text, size = 9),
                        axis.title       = element_text(color = col_text, size = 10),
                        axis.title.y     = element_text(margin = margin(r = 10)),
                        axis.title.x     = element_text(margin = margin(t = 8)),
                        plot.margin      = margin(15, 15, 10, 10)
                )
        
        return(p)
}


# PREDICTOR_SPEC ----------------------------------------------------------


# Full predictor specification with metadata
PREDICTOR_SPEC <- list(
        # Typology characteristics
        n_types = list(
                label = "Number of Types",
                type = "numeric",
                is_integer = TRUE,
                min = 2, max = 50, value = 10,
                group = "Typology"
        ),
        variables = list(
                label = "Number of Variables",
                type = "numeric",
                is_integer = TRUE,
                min = 1, max = 100, value = 5,
                group = "Typology"
        ),
        env_asw = list(
                label = "Environmental ASW",
                type = "numeric",
                min = -1, max = 1, value = 0.3, step = 0.01,
                group = "Typology"
        ),
        fuzzy_npe = list(
                label = "Fuzzy NPE",
                type = "numeric",
                min = 0, max = 1, value = 0.5, step = 0.01,
                group = "Typology",
                conditional = c("fuzzy_mantel")
        ),
        
        # Spatial characteristics
        samples = list(
                label = "Number of Samples",
                type = "numeric",
                is_integer = TRUE,
                min = 1, max = 10000, value = 500,
                group = "Spatial"
        ),
        min_distance = list(
                label = "Min Distance (km)",
                type = "numeric",
                min = 1, max = 10000, value = 10, step = 1,
                group = "Spatial"
        ),
        mean_distance = list(
                label = "Mean Distance (km)",
                type = "numeric",
                min = 1, max = 5000, value = 500, step = 10,
                group = "Spatial"
        ),
        median_distance = list(
                label = "Median Distance (km)",
                type = "numeric",
                min = 1, max = 5000, value = 400, step = 10,
                group = "Spatial"
        ),
        max_distance = list(
                label = "Max Distance (km)",
                type = "numeric",
                min = 1, max = 10000, value = 2000, step = 10,
                group = "Spatial"
        ),
        median_latitude = list(
                label = "Median Latitude",
                type = "numeric",
                min = -90, max = 90, value = 50, step = 0.1,
                group = "Spatial"
        ),
        median_longitude = list(
                label = "Median Longitude",
                type = "numeric",
                min = -180, max = 180, value = 10, step = 0.1,
                group = "Spatial"
        ),
        max_latitude = list(
                label = "Max Latitude",
                type = "numeric",
                min = -90, max = 90, value = 70, step = 0.1,
                group = "Spatial"
        ),
        max_longitude = list(
                label = "Max Longitude",
                type = "numeric",
                min = -180, max = 180, value = 30, step = 0.1,
                group = "Spatial"
        ),
        min_latitude = list(
                label = "Min Latitude",
                type = "numeric",
                min = -90, max = 90, value = 35, step = 0.1,
                group = "Spatial"
        ),
        min_longitude = list(
                label = "Min Longitude",
                type = "numeric",
                min = -180, max = 180, value = -10, step = 0.1,
                group = "Spatial"
        ),
        
        # Variation partitioning
        bio = list(
                label = "Biotic Variation (Fraction)",
                type = "numeric",
                min = 0, max = 1, value = .05, step = .1,
                group = "Variation"
        ),
        env = list(
                label = "Environmental Variation (Fraction)",
                type = "numeric",
                min = 0, max = 1, value = .1, step = .1,
                group = "Variation"
        ),
        space = list(
                label = "Spatial Variation (Fraction)",
                type = "numeric",
                min = 0, max = 1, value = .05, step = .1,
                group = "Variation"
        ),
        stochastic = list(
                label = "Stochastic Variation (Fraction)",
                type = "numeric",
                min = 0, max = 1, value = .8, step = .1,
                group = "Variation"
        ),
        
        # Taxonomic resolution
        n_taxa = list(
                label = "Number of Taxa",
                type = "numeric",
                is_integer = TRUE,
                min = 1, max = 10000, value = 25,
                group = "Taxonomic"
        ),
        species_rank = list(
                label = "Species Rank (Fraction)",
                type = "numeric",
                min = 0, max = 1, value = .7, step = .1,
                group = "Taxonomic"
        ),
        genus_rank = list(
                label = "Genus Rank (%)",
                type = "numeric",
                min = 0, max = 1, value = 0.2, step = .1,
                group = "Taxonomic"
        ),
        family_rank = list(
                label = "Family Rank (%)",
                type = "numeric",
                min = 0, max = 1, value = .08, step = .1,
                group = "Taxonomic"
        ),
        higher_rank = list(
                label = "Higher Rank (%)",
                type = "numeric",
                min = 0, max = 1, value = .02, step = .1,
                group = "Taxonomic"
        )
)

# Group predictors
PREDICTOR_GROUPS <- unique(sapply(PREDICTOR_SPEC, function(x) x$group))

# Custom Theme — Freshwater Ecology Aesthetic
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
  ")


# UI ----------------------------------------------------------------------
ui <- page_navbar(
        title = tags$span(
                tags$span(class = "status-indicator status-ready"),
                "QRF Predictor"
        ),
        theme = custom_theme,
        fillable = TRUE,
        
        # Main prediction tab
        nav_panel(
                title = "Predict",
                icon = icon("chart-line"),
                
                layout_sidebar(
                        sidebar = sidebar(
                                title = "Parameters",
                                width = 380,
                                
                                # Model selection
                                card(
                                        card_header("Model Selection"),
                                        
                                        # Taxonomic group selection
                                        radioButtons(
                                                "taxon_group",
                                                "Taxonomic Group",
                                                choices = setNames(
                                                        names(TAXONOMIC_GROUPS),
                                                        sapply(TAXONOMIC_GROUPS, function(x) x$label)
                                                ),
                                                selected = "macroinvertebrates"
                                        ),
                                        
                                        hr(),
                                        
                                        selectInput(
                                                "metric",
                                                "Evaluation Metric",
                                                choices = NULL,
                                                selected = NULL
                                        ),
                                        
                                        checkboxGroupInput(
                                                "quantiles",
                                                "Prediction Quantiles",
                                                choices = c(
                                                        "5%" = 0.05,
                                                        "10%" = 0.10,
                                                        "25%" = 0.25,
                                                        "50% (Median)" = 0.50,
                                                        "75%" = 0.75,
                                                        "90%" = 0.90,
                                                        "95%" = 0.95
                                                ),
                                                selected = c(0.10, 0.50, 0.90),
                                                inline = TRUE
                                        ),
                                        
                                        hr(),
                                        
                                        # Observed value input
                                        div(
                                                class = "observed-value-section",
                                                numericInput(
                                                        "observed_value",
                                                        tags$span(
                                                                icon("crosshairs"),
                                                                " Your Observed Value"
                                                        ),
                                                        value = NA,
                                                        step = 0.001
                                                ),
                                                tags$small(
                                                        class = "text-muted",
                                                        "Optional. Enter your typology's actual metric ",
                                                        "value to see where it falls in the predicted ",
                                                        "distribution."
                                                )
                                        ),
                                        
                                        hr(),
                                        
                                        # Imputation method for unspecified parameters
                                        radioButtons(
                                                "impute_method",
                                                "Unspecified Parameters",
                                                choices = c(
                                                        "Hard-coded defaults" = "defaults",
                                                        "Training medians" = "marginal",
                                                        "Conditional averages (KNN)" = "conditional"
                                                ),
                                                selected = "defaults"
                                        ),
                                        
                                        # KNN settings (only shown when conditional is selected)
                                        conditionalPanel(
                                                condition = "input.impute_method == 'conditional'",
                                                sliderInput(
                                                        "knn_k",
                                                        "Neighbors (k)",
                                                        min = 5, max = 200, value = 50, step = 5
                                                ),
                                                uiOutput("training_data_status")
                                        )
                                ),
                                
                                # Parameter toggles
                                card(
                                        card_header("Active Parameters"),
                                        p(class = "text-muted small", 
                                          "Toggle which parameters to include in prediction. Unselected parameters use the method chosen above."),
                                        
                                        navset_pill(
                                                id = "param_tabs",
                                                
                                                nav_panel(
                                                        "Typology",
                                                        div(class = "mt-3",
                                                            uiOutput("typology_toggles")
                                                        )
                                                ),
                                                nav_panel(
                                                        "Spatial",
                                                        div(class = "mt-3",
                                                            uiOutput("spatial_toggles")
                                                        )
                                                ),
                                                nav_panel(
                                                        "Variation Partitioning",
                                                        div(class = "mt-3",
                                                            uiOutput("variation_toggles")
                                                        )
                                                ),
                                                nav_panel(
                                                        "Taxonomic",
                                                        div(class = "mt-3",
                                                            uiOutput("taxonomic_toggles")
                                                        )
                                                )
                                        )
                                ),
                                
                                # Action buttons
                                div(
                                        class = "d-grid gap-2 mt-3",
                                        actionButton("predict", "Generate Prediction", 
                                                     class = "btn-primary btn-lg",
                                                     icon = icon("wand-magic-sparkles")),
                                        actionButton("reset", "Reset to Defaults", 
                                                     class = "btn-outline-secondary",
                                                     icon = icon("rotate-left"))
                                ),
                                # Remove cache 
                                card(
                                        card_header("Cache"),
                                        uiOutput("cache_status"),
                                        div(
                                                class = "d-grid gap-2 mt-2",
                                                actionButton("clear_cache", "Clear Cached Files",
                                                             class = "btn-outline-secondary",
                                                             icon = icon("trash-can"))
                                        ),
                                        tags$small(
                                                class = "text-muted",
                                                "Downloaded models and training data are cached locally for faster ",
                                                "access. Clearing the cache forces fresh downloads on next use."
                                        )
                                )
                        ),
                        
                        # Main content area
                        layout_columns(
                                col_widths = c(12, 8, 4),
                                
                                # Prediction result card (now includes density plot)
                                card(
                                        card_header(
                                                class = "d-flex justify-content-between align-items-center",
                                                "Prediction Result",
                                                uiOutput("model_status")
                                        ),
                                        
                                        uiOutput("prediction_display"),
                                        
                                        # Density plot section
                                        div(
                                                class = "density-plot-container",
                                                div(class = "density-plot-header",
                                                    icon("chart-area"),
                                                    " Predicted Distribution"
                                                ),
                                                div(class = "density-plot-subtitle",
                                                    "Density estimated from a dense quantile grid (199 points)"
                                                ),
                                                plotOutput("density_plot", height = "300px"),
                                                uiOutput("observed_percentile_display")
                                        )
                                ),
                                
                                # Parameter inputs
                                card(
                                        card_header("Parameter Values"),
                                        div(
                                                style = "max-height: 500px; overflow-y: auto;",
                                                uiOutput("parameter_inputs")
                                        )
                                ),
                                
                                # Quick stats
                                card(
                                        card_header("Model Info"),
                                        uiOutput("model_info")
                                )
                        )
                )
        ),
        
        # Batch prediction tab
        nav_panel(
                title = "Batch",
                icon = icon("table"),
                
                layout_columns(
                        col_widths = c(4, 8),
                        
                        card(
                                card_header("Upload Data"),
                                fileInput("batch_file", "CSV File", 
                                          accept = c(".csv", ".CSV")),
                                
                                radioButtons(
                                        "batch_taxon_group",
                                        "Taxonomic Group",
                                        choices = setNames(
                                                names(TAXONOMIC_GROUPS),
                                                sapply(TAXONOMIC_GROUPS, function(x) x$label)
                                        ),
                                        selected = "macroinvertebrates",
                                        inline = TRUE
                                ),
                                
                                selectInput(
                                        "batch_metric",
                                        "Evaluation Metric",
                                        choices = NULL,
                                        selected = NULL
                                ),
                                
                                checkboxGroupInput(
                                        "batch_quantiles",
                                        "Quantiles to Predict",
                                        choices = c("10%" = 0.1, "25%" = 0.25, "50%" = 0.5, 
                                                    "75%" = 0.75, "90%" = 0.9),
                                        selected = c(0.1, 0.5, 0.9)
                                ),
                                
                                hr(),
                                
                                radioButtons(
                                        "batch_impute_method",
                                        "Missing Column Handling",
                                        choices = c(
                                                "Hard-coded defaults" = "defaults",
                                                "Training medians" = "marginal",
                                                "Conditional averages (KNN)" = "conditional"
                                        ),
                                        selected = "defaults"
                                ),
                                
                                conditionalPanel(
                                        condition = "input.batch_impute_method == 'conditional'",
                                        sliderInput(
                                                "batch_knn_k",
                                                "Neighbors (k)",
                                                min = 5, max = 200, value = 50, step = 5
                                        )
                                ),
                                
                                uiOutput("batch_column_status"),
                                
                                div(
                                        class = "d-grid gap-2 mt-3",
                                        actionButton("batch_predict", "Run Batch Prediction",
                                                     class = "btn-primary",
                                                     icon = icon("play")),
                                        downloadButton("download_results", "Download Results",
                                                       class = "btn-outline-secondary")
                                )
                        ),
                        
                        card(
                                card_header("Results"),
                                DTOutput("batch_results")
                        )
                )
        ),
        
        # About tab
        nav_panel(
                title = "About",
                icon = icon("circle-info"),
                
                layout_columns(
                        col_widths = c(6, 6),
                        
                        card(
                                card_header("About This App"),
                                markdown("
This application provides predictions from **Quantile Random Forest** models 
fitted on freshwater typology evaluation metrics.

### Taxonomic Groups

Models are available for four taxonomic groups:

- 🌿 **Diatoms** — Microscopic algae
- 🐟 **Fishes** — Fish assemblages
- 🦗 **Macroinvertebrates** — Benthic invertebrates
- 🌱 **Macrophytes** — Aquatic plants

### Features

- **Single predictions** with uncertainty quantiles
- **Predicted distribution visualization** with observed value marking
- **Batch predictions** from CSV uploads
- **Flexible parameter selection** — use only the parameters you have
- **Multiple metrics** supported per taxonomic group
- **Conditional averaging** — unspecified parameters are imputed using 
  KNN on training data, conditioned on your specified values

### How It Works

1. Select a taxonomic group
2. Select an evaluation metric
3. Choose how unspecified parameters should be filled:
   - *Hard-coded defaults*: fixed values from the predictor specification
   - *Training medians*: marginal medians from the training dataset
   - *Conditional averages (KNN)*: averages from the k nearest neighbors
     in the training data, based on your specified parameter values
4. Toggle which parameters you want to specify
5. Enter values for active parameters
6. Optionally enter your observed metric value
7. Click predict to get quantile predictions and a density visualization

### Distribution Visualization

After prediction, the app displays an estimated density curve by:

1. Predicting on a dense quantile grid (199 points from Q0.5 to Q99.5)
2. Numerically differentiating the quantile function: density ≈ Δp / Δx
3. Applying loess smoothing for visual clarity
4. Marking the user-selected quantiles and (optionally) the observed value
5. Reporting the percentile rank of the observed value in the predicted
   distribution

### Unspecified Parameter Imputation

When using **Conditional Averages**, the app:

1. Takes the variables you've specified and their values
2. Standardises these variables using the training data (z-scores)
3. Computes Euclidean distance to every row in the training data
4. Selects the *k* nearest neighbors
5. Averages the unspecified variables across those neighbors

This means the imputed values reflect realistic covariance structure 
in the training data rather than marginal or hard-coded defaults.

### Directory Structure

Models and training data are fetched on demand from Zenodo and cached locally. 
No local file setup is required. The model record is available at https://doi.org/10.5281/zenodo.19554225.

Training data files (`training_data.rds`) should be a `data.frame` 
containing the predictor columns used by the models. Alternatively, 
metric-specific training data (`<metric>_training_data.rds`) can be 
provided and will take priority.
        ")
                        ),
                        
                        card(
                                card_header("Parameter Reference"),
                                div(
                                        style = "max-height: 400px; overflow-y: auto;",
                                        DTOutput("param_reference")
                                )
                        )
                )
        )
)

# Server ------------------------------------------------------------------

server <- function(input, output, session) {
        
        # Reactive: Store loaded model
        loaded_model <- reactiveVal(NULL)
        model_predictors <- reactiveVal(NULL)
        training_medians <- reactiveVal(NULL)
        
        # Reactive: Store training data
        loaded_training_data <- reactiveVal(NULL)
        
        # Update available metrics when taxonomic group changes
        observeEvent(input$taxon_group, {
                req(input$taxon_group)
                cat("Taxonomic group changed to:", input$taxon_group, "\n")
                metrics <- get_available_metrics(input$taxon_group)
                cat("Available metrics:", paste(metrics, collapse = ", "), "\n")
                updateSelectInput(session, "metric", 
                                  choices = metrics,
                                  selected = if(length(metrics) > 0) metrics[1] else NULL)
        }, ignoreNULL = TRUE)
        
        # Update batch metrics when batch taxonomic group changes
        observeEvent(input$batch_taxon_group, {
                req(input$batch_taxon_group)
                metrics <- get_available_metrics(input$batch_taxon_group)
                updateSelectInput(session, "batch_metric",
                                  choices = metrics,
                                  selected = if(length(metrics) > 0) metrics[1] else NULL)
        }, ignoreNULL = TRUE)
        
        # Reactive trigger so the status updates after a clear
        cache_invalidator <- reactiveVal(0)
        
        output$cache_status <- renderUI({
                cache_invalidator()  # take a dependency so this re-runs after clearing
                info <- cache_summary()
                if (info$n_files == 0) {
                        div(
                                class = "small",
                                style = "color: var(--text-light);",
                                icon("folder-open"),
                                " Cache is empty"
                        )
                } else {
                        div(
                                class = "small",
                                style = "color: var(--accent);",
                                icon("database"),
                                sprintf(" %d file%s cached (%.1f MB)",
                                        info$n_files,
                                        if (info$n_files == 1) "" else "s",
                                        info$total_mb)
                        )
                }
        })

        # Remove cache  -----------------------------------------------------------
        observeEvent(input$clear_cache, {
                # Show a confirmation modal rather than wiping immediately
                showModal(modalDialog(
                        title = "Clear cached files?",
                        paste0("This will delete all downloaded models and training data. ",
                               "The next prediction will re-download them from Zenodo, ",
                               "which may take several minutes."),
                        footer = tagList(
                                modalButton("Cancel"),
                                actionButton("confirm_clear_cache", "Clear cache",
                                             class = "btn-danger")
                        )
                ))
        })
        
        observeEvent(input$confirm_clear_cache, {
                removeModal()
                result <- clear_cache()
                
                # Reset in-memory state so the app doesn't hold stale references
                loaded_model(NULL)
                model_predictors(NULL)
                loaded_training_data(NULL)
                prediction_result(NULL)
                
                # Trigger the cache status UI to refresh
                cache_invalidator(cache_invalidator() + 1)
                
                if (result$n_files == 0) {
                        showNotification("Cache was already empty.", type = "message")
                } else {
                        showNotification(
                                sprintf("Cleared %d file%s (%.1f MB)",
                                        result$n_files,
                                        if (result$n_files == 1) "" else "s",
                                        result$total_mb),
                                type = "message",
                                duration = 4
                        )
                }
        })
        
        
        # Load model AND training data when metric OR taxonomic group changes
        observeEvent(c(input$metric, input$taxon_group), {
                req(input$metric, input$taxon_group)
                req(nchar(input$metric) > 0, nchar(input$taxon_group) > 0)
                
                # ── Fetch model file from Zenodo (or cache) ──
                model_filename <- build_remote_filename(
                        input$taxon_group, input$metric, type = "model"
                )
                
                withProgress(message = "Loading model", value = 0.3, {
                        incProgress(0.2, detail = "Fetching from Zenodo...")
                        model_file <- fetch_remote_file(model_filename)
                        
                        if (is.null(model_file)) {
                                loaded_model(NULL)
                                model_predictors(NULL)
                                showNotification(
                                        paste("Could not fetch model:", model_filename),
                                        type = "error", duration = 5
                                )
                                return()
                        }
                        incProgress(0.3, detail = "Unbundling...")
                        tryCatch({
                                bundled <- readRDS(model_file)
                                model <- unbundle(bundled)
                                loaded_model(model)
                                
                                recipe <- extract_recipe(model)
                                if (is.null(recipe)) stop("Could not extract recipe from model")
                                
                                predictors <- tryCatch({
                                        recipe$var_info$variable[recipe$var_info$role == "predictor"]
                                }, error = function(e) {
                                        names(recipe$template)[names(recipe$template) != "value"]
                                })
                                
                                model_predictors(predictors)
                                training_medians(sapply(PREDICTOR_SPEC, function(x) x$value))
                                
                                showNotification(paste("Model loaded:", input$metric),
                                                 type = "message", duration = 3)
                        }, error = function(e) {
                                loaded_model(NULL)
                                model_predictors(NULL)
                                showNotification(paste("Error loading model:", e$message), type = "error")
                        })
                        
                        # ── Fetch training data ──
                        incProgress(0.2, detail = "Fetching training data...")
                        td_filename <- build_remote_filename(
                                input$taxon_group, type = "training_data"
                        )
                        td_file <- fetch_remote_file(td_filename)
                        if (!is.null(td_file)) {
                                tryCatch({
                                        td <- readRDS(td_file)
                                        if (inherits(td, "data.table")) 
                                                td <- as.data.frame(td)
                                        loaded_training_data(td)
                                        preds <- model_predictors()
                                        if (!is.null(preds)) {
                                                available_cols <- intersect(preds, names(td))
                                                real_medians <- sapply(available_cols, function(v) median(td[[v]], na.rm = TRUE))
                                                all_medians <- sapply(PREDICTOR_SPEC, function(x) x$value)
                                                all_medians[names(real_medians)] <- real_medians
                                                training_medians(all_medians)
                                        }
                                }, error = function(e) {
                                        loaded_training_data(NULL)
                                })
                        } else {
                                loaded_training_data(NULL)
                        }
                        # } else {
                        #         loaded_training_data(NULL)
                        # }
                        
                        incProgress(0.3, detail = "Done")
                })
        }, ignoreInit = TRUE)

        # Training data status indicator
        output$training_data_status <- renderUI({
                td <- loaded_training_data()
                if (!is.null(td)) {
                        div(
                                class = "small mt-2",
                                style = "color: var(--accent);",
                                icon("check-circle"),
                                sprintf(" Training data loaded (%s rows)", 
                                        format(nrow(td), big.mark = ","))
                        )
                } else {
                        div(
                                class = "small mt-2",
                                style = "color: var(--danger);",
                                icon("triangle-exclamation"),
                                " No training data found. Will fall back to defaults.",
                                br(),
                                tags$span(
                                        class = "text-muted",
                                        style = "font-size: 0.78rem;",
                                        "Training data could not be fetched from Zenodo. Check your internet connection."
                                )
                        )
                }
        })
        
        # Model status indicator
        output$model_status <- renderUI({
                if (!is.null(loaded_model())) {
                        taxon_label <- TAXONOMIC_GROUPS[[input$taxon_group]]$label
                        tags$span(
                                class = "badge bg-success",
                                icon("check"), paste(taxon_label, "- Ready")
                        )
                } else {
                        tags$span(
                                class = "badge bg-warning",
                                icon("exclamation-triangle"), "No Model"
                        )
                }
        })
        
        # Generate toggle checkboxes for each group
        create_toggle_ui <- function(group_name) {
                params <- names(PREDICTOR_SPEC)[sapply(PREDICTOR_SPEC, function(x) x$group == group_name)]
                
                tagList(
                        lapply(params, function(p) {
                                spec <- PREDICTOR_SPEC[[p]]
                                
                                is_available <- TRUE
                                if (!is.null(spec$conditional)) {
                                        is_available <- input$metric %in% spec$conditional
                                }
                                
                                if (is_available) {
                                        checkboxInput(
                                                inputId = paste0("toggle_", p),
                                                label = spec$label,
                                                value = TRUE
                                        )
                                }
                        })
                )
        }
        
        output$typology_toggles <- renderUI({ create_toggle_ui("Typology") })
        output$spatial_toggles <- renderUI({ create_toggle_ui("Spatial") })
        output$variation_toggles <- renderUI({ create_toggle_ui("Variation") })
        output$taxonomic_toggles <- renderUI({ create_toggle_ui("Taxonomic") })
        
        # Generate parameter input fields
        output$parameter_inputs <- renderUI({
                required_preds <- model_predictors()
                
                groups <- lapply(PREDICTOR_GROUPS, function(group) {
                        params <- names(PREDICTOR_SPEC)[sapply(PREDICTOR_SPEC, function(x) x$group == group)]
                        
                        params <- params[sapply(params, function(p) {
                                spec <- PREDICTOR_SPEC[[p]]
                                
                                toggle_id <- paste0("toggle_", p)
                                is_active <- tryCatch(
                                        isTRUE(input[[toggle_id]]) || is.null(input[[toggle_id]]),
                                        error = function(e) TRUE
                                )
                                
                                is_available <- TRUE
                                if (!is.null(spec$conditional)) {
                                        is_available <- isTRUE(input$metric %in% spec$conditional)
                                }
                                in_model <- is.null(required_preds) || p %in% required_preds
                                
                                is_active && is_available && in_model
                        })]
                        
                        if (length(params) == 0) return(NULL)
                        
                        div(
                                div(class = "parameter-group-title", group),
                                lapply(params, function(p) {
                                        spec <- PREDICTOR_SPEC[[p]]
                                        
                                        if (spec$type == "numeric") {
                                                # Preserve current value if user already set one
                                                current_val <- isolate(input[[paste0("param_", p)]])
                                                use_val <- if (!is.null(current_val)) current_val else spec$value
                                                
                                                numericInput(
                                                        inputId = paste0("param_", p),
                                                        label = spec$label,
                                                        value = use_val,
                                                        min = spec$min,
                                                        max = spec$max,
                                                        step = spec$step %||% 1
                                                )
                                        } else if (spec$type == "select") {
                                                selectInput(
                                                        inputId = paste0("param_", p),
                                                        label = spec$label,
                                                        choices = spec$choices,
                                                        selected = spec$value
                                                )
                                        }
                                })
                        )
                })
                
                tagList(groups)
        })
        
        # Model info display
        output$model_info <- renderUI({
                if (is.null(loaded_model())) {
                        return(div(
                                class = "text-center p-4",
                                style = "color: var(--text-light);",
                                icon("folder-open", class = "fa-3x mb-3"),
                                p("No model loaded"),
                                p(class = "small", "Select a taxonomic group and metric")
                        ))
                }
                
                preds <- model_predictors()
                taxon_info <- TAXONOMIC_GROUPS[[input$taxon_group]]
                td <- loaded_training_data()
                
                div(
                        # Taxonomic group indicator
                        div(
                                class = "text-center mb-3 p-2",
                                style = "background: rgba(46, 107, 138, 0.06); border-radius: 8px;",
                                icon(taxon_info$icon, class = "fa-2x mb-2", 
                                     style = "color: var(--water-blue);"),
                                div(class = "fw-bold", 
                                    style = "color: var(--river-deep);", 
                                    taxon_info$label)
                        ),
                        
                        div(class = "metric-card",
                            div(class = "metric-label", "Predictors"),
                            div(class = "h3", 
                                style = "color: var(--river-deep); font-family: 'Instrument Serif', serif;",
                                length(preds))
                        ),
                        
                        if (!is.null(td)) {
                                div(class = "metric-card",
                                    div(class = "metric-label", "Training Samples"),
                                    div(class = "h4",
                                        style = "color: var(--water-blue); font-family: 'Instrument Serif', serif;",
                                        format(nrow(td), big.mark = ","))
                                )
                        },
                        
                        hr(),
                        div(
                                class = "small",
                                tags$strong(style = "color: var(--river-deep);", "Variables:"),
                                p(class = "text-muted", 
                                  paste(preds, collapse = ", "))
                        )
                )
        })
        
        # Prediction result (stores extra info about imputation + density data)
        prediction_result <- reactiveVal(NULL)
        
        
        # Core prediction logic with imputation method selection
        
        observeEvent(input$predict, {
                req(loaded_model())
                
                model <- loaded_model()
                preds <- model_predictors()
                medians <- training_medians()
                td <- loaded_training_data()
                
                # ── Separate specified vs unspecified predictors ──
                specified_vals <- list()
                unspecified_vars <- c()
                
                for (p in preds) {
                        param_input <- input[[paste0("param_", p)]]
                        
                        toggle_input <- tryCatch(
                                isTRUE(input[[paste0("toggle_", p)]]),
                                error = function(e) TRUE
                        )
                        if (is.null(input[[paste0("toggle_", p)]])) {
                                toggle_input <- TRUE
                        }
                        
                        if (!is.null(param_input) && toggle_input) {
                                specified_vals[[p]] <- param_input
                        } else {
                                unspecified_vars <- c(unspecified_vars, p)
                        }
                }
                
                # ── Resolve unspecified variable values ──
                impute_method <- input$impute_method
                imputed_vals <- list()
                
                if (length(unspecified_vars) > 0) {
                        if (impute_method == "conditional" && !is.null(td)) {
                                # KNN conditional averaging
                                k <- input$knn_k %||% 50
                                cond_avg <- compute_conditional_averages(
                                        training_data = td,
                                        specified_vals = specified_vals,
                                        unspecified_vars = unspecified_vars,
                                        k = k
                                )
                                for (v in unspecified_vars) {
                                        imputed_vals[[v]] <- if (v %in% names(cond_avg) && !is.na(cond_avg[[v]])) {
                                                cond_avg[[v]]
                                        } else {
                                                # Fallback
                                                if (!is.null(medians) && v %in% names(medians)) medians[[v]]
                                                else if (v %in% names(PREDICTOR_SPEC)) PREDICTOR_SPEC[[v]]$value
                                                else 0
                                        }
                                }
                                
                        } else if (impute_method == "marginal" && !is.null(td)) {
                                # Marginal medians from training data
                                for (v in unspecified_vars) {
                                        imputed_vals[[v]] <- if (v %in% names(td)) {
                                                median(td[[v]], na.rm = TRUE)
                                        } else if (!is.null(medians) && v %in% names(medians)) {
                                                medians[[v]]
                                        } else if (v %in% names(PREDICTOR_SPEC)) {
                                                PREDICTOR_SPEC[[v]]$value
                                        } else {
                                                0
                                        }
                                }
                                
                        } else {
                                # Hard-coded defaults
                                for (v in unspecified_vars) {
                                        imputed_vals[[v]] <- if (!is.null(medians) && v %in% names(medians)) {
                                                medians[[v]]
                                        } else if (v %in% names(PREDICTOR_SPEC)) {
                                                PREDICTOR_SPEC[[v]]$value
                                        } else {
                                                0
                                        }
                                }
                        }
                }
                
                # ── Compositional constraint enforcement ──
                for (grp in COMPOSITIONAL_GROUPS) {
                        spec_in_grp <- intersect(grp, names(specified_vals))
                        unspec_in_grp <- intersect(grp, unspecified_vars)
                        
                        if (length(unspec_in_grp) > 0) {
                                sum_specified <- sum(unlist(specified_vals[spec_in_grp]))
                                remaining <- max(0, 1 - sum_specified)
                                raw_imputed <- unlist(imputed_vals[unspec_in_grp])
                                sum_raw <- sum(raw_imputed)
                                
                                if (sum_raw > 0) {
                                        scaled <- (raw_imputed / sum_raw) * remaining
                                } else {
                                        scaled <- rep(remaining / length(unspec_in_grp), length(unspec_in_grp))
                                }
                                
                                for (v in unspec_in_grp) {
                                        imputed_vals[[v]] <- scaled[[v]]
                                }
                        }
                }
                
                # ── Enforce ordered spatial constraints ──
                for (oc in ORDERED_CONSTRAINTS) {
                        chain    <- oc$chain
                        interior <- oc$interior
                        endpoints <- oc$endpoints
                        
                        # Gather current values
                        get_val <- function(v) {
                                if (v %in% names(specified_vals)) specified_vals[[v]]
                                else if (v %in% names(imputed_vals)) imputed_vals[[v]]
                                else NA
                        }
                        
                        lo_var <- endpoints[1]
                        hi_var <- endpoints[2]
                        lo_val <- get_val(lo_var)
                        hi_val <- get_val(hi_var)
                        
                        # Identify which interior vars are imputed vs specified
                        imputed_interior  <- intersect(interior, unspecified_vars)
                        specified_interior <- setdiff(interior, unspecified_vars)
                        
                        # Case 1: All interior vars are imputed and both endpoints exist
                        #         → place imputed interior at the midpoint
                        if (length(imputed_interior) == length(interior) &&
                            !is.na(lo_val) && !is.na(hi_val)) {
                                midpoint <- (lo_val + hi_val) / 2
                                for (v in imputed_interior) {
                                        imputed_vals[[v]] <- midpoint
                                }
                                
                                # Case 2: Some interior specified, some imputed
                                #         → set imputed ones to match the specified one
                        } else if (length(imputed_interior) > 0 && length(specified_interior) > 0) {
                                ref_val <- get_val(specified_interior[1])
                                for (v in imputed_interior) {
                                        imputed_vals[[v]] <- ref_val
                                }
                        }
                        
                        # Safety-net clamping pass (only adjusts imputed values)
                        chain_vals <- sapply(chain, get_val)
                        
                        for (i in 2:length(chain)) {
                                if (!is.na(chain_vals[i]) && !is.na(chain_vals[i - 1])) {
                                        if (chain_vals[i] < chain_vals[i - 1] && chain[i] %in% unspecified_vars) {
                                                chain_vals[i] <- chain_vals[i - 1]
                                                imputed_vals[[chain[i]]] <- chain_vals[i]
                                        }
                                }
                        }
                        for (i in (length(chain) - 1):1) {
                                if (!is.na(chain_vals[i]) && !is.na(chain_vals[i + 1])) {
                                        if (chain_vals[i] > chain_vals[i + 1] && chain[i] %in% unspecified_vars) {
                                                chain_vals[i] <- chain_vals[i + 1]
                                                imputed_vals[[chain[i]]] <- chain_vals[i]
                                        }
                                }
                        }
                }
                
                # ── Build prediction data.frame ──
                new_data <- data.frame(row.names = 1)
                for (p in preds) {
                        if (p %in% names(specified_vals)) {
                                val <- specified_vals[[p]]
                        } else {
                                val <- imputed_vals[[p]]
                        }
                        
                        if (isTRUE(PREDICTOR_SPEC[[p]]$is_integer)) {
                                val <- round(val)
                                if (!(p %in% names(specified_vals))) {
                                        imputed_vals[[p]] <- val
                                }
                        }
                        new_data[[p]] <- val
                }
                
                # ── Run prediction ──
                tryCatch({
                        ranger_model <- extract_fit_engine(model)
                        recipe <- extract_recipe(model)
                        
                        recipe_vars <- recipe$var_info$variable
                        missing_vars <- setdiff(recipe_vars, names(new_data))
                        if (length(missing_vars) > 0) {
                                for (mv in missing_vars) {
                                        new_data[[mv]] <- NA
                                }
                        }
                        
                        prepped_data <- bake(recipe, new_data = as.data.table(new_data))
                        
                        # User-selected quantiles for the badge display
                        quantiles <- as.numeric(input$quantiles)
                        if (length(quantiles) == 0) {
                                quantiles <- c(0.1, 0.5, 0.9)
                        }
                        quantiles <- sort(quantiles)
                        
                        qpred <- predict(
                                ranger_model,
                                data = prepped_data,
                                type = "quantiles",
                                quantiles = quantiles
                        )
                        
                        # ── Dense quantile grid for density estimation ──
                        dens_df <- tryCatch({
                                # estimate_density_from_qrf(
                                #         ranger_model = ranger_model,
                                #         prepped_data = prepped_data,
                                #         from = 0.005,
                                #         to = 0.995,
                                #         n_quantiles = 199,
                                #         smooth = TRUE,
                                #         smooth_span = 0.35
                                # )
                                estimate_density_from_qrf(ranger_model = ranger_model,
                                                          prepped_data = prepped_data,
                                                          n_samples = 500,
                                                          bw_method = "nrd0")
                        }, error = function(e) {
                                cat("  Density estimation failed:", e$message, "\n")
                                NULL
                        })
                        
                        # ── Compute observed percentile if provided ──
                        obs_val <- input$observed_value
                        obs_percentile <- NULL
                        
                        if (!is.null(obs_val) && !is.na(obs_val) && is.finite(obs_val)) {
                                # Predict at a very fine grid to find the percentile
                                obs_quantiles <- seq(0.001, 0.999, by = 0.001)
                                obs_qpred <- predict(
                                        ranger_model,
                                        data = prepped_data,
                                        type = "quantiles",
                                        quantiles = obs_quantiles
                                )
                                obs_qvals <- as.numeric(obs_qpred$predictions[1, ])
                                
                                # Find the percentile: fraction of predicted quantile 
                                # values that are <= the observed value
                                obs_percentile <- tryCatch({
                                        approx(obs_qvals, obs_quantiles, 
                                               xout = obs_val, rule = 2)$y
                                }, error = function(e) {
                                        # Fallback: simple empirical CDF
                                        mean(obs_qvals <= obs_val)
                                })
                        }
                        
                        prediction_result(list(
                                quantiles = quantiles,
                                values = as.numeric(qpred$predictions[1, ]),
                                input_data = new_data,
                                impute_method = impute_method,
                                specified_vars = names(specified_vals),
                                imputed_vars = unspecified_vars,
                                imputed_vals = imputed_vals,
                                density_data = dens_df,
                                observed_value = if (!is.null(obs_val) && !is.na(obs_val)) obs_val else NULL,
                                observed_percentile = obs_percentile
                        ))
                        
                }, error = function(e) {
                        showNotification(
                                paste("Prediction error:", e$message),
                                type = "error"
                        )
                        prediction_result(NULL)
                })
        })
        
        
        # Display prediction (text + badges)
        
        output$prediction_display <- renderUI({
                res <- prediction_result()
                
                if (is.null(res)) {
                        return(div(
                                class = "text-center p-5",
                                icon("chart-line", class = "fa-4x mb-3", 
                                     style = "color: var(--border);"),
                                h4(style = "color: var(--text-light); font-family: 'Instrument Serif', serif;", 
                                   "No prediction yet"),
                                p(style = "color: var(--text-light);", 
                                  "Configure parameters and click 'Generate Prediction'")
                        ))
                }
                
                median_idx <- which.min(abs(res$quantiles - 0.5))
                median_val <- res$values[median_idx]
                
                # Build imputation summary
                impute_summary <- NULL
                if (length(res$imputed_vars) > 0) {
                        method_label <- switch(
                                res$impute_method,
                                "conditional" = "Conditional averages (KNN)",
                                "marginal" = "Training medians",
                                "defaults" = "Hard-coded defaults"
                        )
                        
                        impute_summary <- div(
                                class = "imputed-summary",
                                div(class = "imputed-header",
                                    icon("sliders"),
                                    sprintf(" %d imputed via %s", 
                                            length(res$imputed_vars), method_label)
                                ),
                                div(
                                        lapply(res$imputed_vars, function(v) {
                                                val <- res$imputed_vals[[v]]
                                                label <- if (v %in% names(PREDICTOR_SPEC)) {
                                                        PREDICTOR_SPEC[[v]]$label
                                                } else {
                                                        v
                                                }
                                                span(
                                                        class = "imputed-var",
                                                        sprintf("%s = %.2f", label, val)
                                                )
                                        })
                                )
                        )
                }
                
                div(
                        class = "text-center py-4",
                        
                        div(
                                class = "prediction-value",
                                sprintf("%.3f", median_val)
                        ),
                        
                        p(style = "color: var(--text-light);", class = "mb-4", "Median Prediction"),
                        
                        div(
                                class = "mb-4",
                                lapply(seq_along(res$quantiles), function(i) {
                                        q <- res$quantiles[i]
                                        v <- res$values[i]
                                        
                                        badge_class <- if (q < 0.4) {
                                                "quantile-badge quantile-low"
                                        } else if (q > 0.6) {
                                                "quantile-badge quantile-high"
                                        } else {
                                                "quantile-badge quantile-mid"
                                        }
                                        
                                        span(
                                                class = badge_class,
                                                sprintf("Q%.0f: %.3f", q * 100, v)
                                        )
                                })
                        ),
                        
                        if (length(res$values) >= 2) {
                                div(
                                        class = "small",
                                        style = "color: var(--text-light);",
                                        sprintf(
                                                "%.0f%% Prediction Interval: [%.3f, %.3f]",
                                                (max(res$quantiles) - min(res$quantiles)) * 100,
                                                min(res$values),
                                                max(res$values)
                                        )
                                )
                        },
                        
                        # Show imputation summary below prediction
                        impute_summary
                )
        })
        
        
        # Density plot
        
        output$density_plot <- renderPlot({
                res <- prediction_result()
                
                if (is.null(res) || is.null(res$density_data)) {
                        # Empty placeholder
                        return(
                                ggplot() +
                                        annotate("text", x = 0.5, y = 0.5,
                                                 label = "Generate a prediction to see the distribution",
                                                 size = 4.5, color = "#7a7a7a",
                                                 family = "sans") +
                                        theme_void() +
                                        theme(plot.background = element_rect(
                                                fill = "transparent", color = NA))
                        )
                }
                
                # Build the metric label from the current selection
                metric_label <- gsub("_", " ", input$metric)
                metric_label <- paste0(
                        toupper(substring(metric_label, 1, 1)),
                        substring(metric_label, 2)
                )
                
                build_density_plot(
                        dens_df = res$density_data,
                        highlight_quantiles = res$quantiles,
                        highlight_values = res$values,
                        observed_value = res$observed_value,
                        metric_label = metric_label
                )
                
        }, bg = "transparent", res = 96)
        
        
        # Observed value percentile display
        
        output$observed_percentile_display <- renderUI({
                res <- prediction_result()
                
                if (is.null(res) || is.null(res$observed_percentile)) {
                        return(NULL)
                }
                
                pct <- res$observed_percentile * 100
                
                # Colour coding based on percentile
                pct_color <- if (pct < 10 || pct > 90) {
                        "#b85c3a"   # extremes — warning
                } else if (pct < 25 || pct > 75) {
                        "#c4a882"   # moderate
                } else {
                        "#1e7a5a"   # central — good
                }
                
                div(
                        class = "text-center mt-2 mb-3",
                        div(
                                class = "percentile-badge",
                                style = paste0(
                                        "background: ", pct_color, "15;",
                                        " color: ", pct_color, ";"
                                ),
                                icon("crosshairs"),
                                sprintf(" Observed value (%.3f) is at the %.1f%s percentile",
                                        res$observed_value, pct,
                                        if (pct == 1) "st" else if (pct == 2) "nd"
                                        else if (pct == 3) "rd" else "th")
                        ),
                        div(
                                class = "small mt-1",
                                style = paste0("color: ", pct_color, ";"),
                                if (pct < 10) {
                                        "Below the 10th percentile — substantially lower than expected"
                                } else if (pct < 25) {
                                        "Below the 25th percentile — lower than typical"
                                } else if (pct > 90) {
                                        "Above the 90th percentile — substantially higher than expected"
                                } else if (pct > 75) {
                                        "Above the 75th percentile — higher than typical"
                                } else {
                                        "Within the interquartile range — consistent with expectations"
                                }
                        )
                )
        })
        
        
        # Reset button
        
        observeEvent(input$reset, {
                for (p in names(PREDICTOR_SPEC)) {
                        spec <- PREDICTOR_SPEC[[p]]
                        updateNumericInput(session, paste0("param_", p), value = spec$value)
                }
                
                updateNumericInput(session, "observed_value", value = NA)
                prediction_result(NULL)
                showNotification("Parameters reset to defaults", type = "message")
        })
        
        
        # Batch prediction --------------------------------------------------------
        
        batch_results <- reactiveVal(NULL)
        
        
        # Batch: column status display (shows after file upload)
        
        output$batch_column_status <- renderUI({
                req(input$batch_file)
                req(input$batch_metric)
                req(input$batch_taxon_group)
                
                # Read just the header to check columns
                tryCatch({
                        data <- fread(input$batch_file$datapath, nrows = 0)
                        csv_cols <- names(data)
                        
                        # Load model to get required predictors
                        model_filename <- build_remote_filename(
                                input$batch_taxon_group, input$batch_metric, type = "model"
                        )
                        model_file <- fetch_remote_file(model_filename)
                        
                        if (is.null(model_file)) return(NULL)
                        
                        bundled <- readRDS(model_file)
                        model <- unbundle(bundled)
                        recipe <- extract_recipe(model)
                        required_preds <- tryCatch({
                                recipe$var_info$variable[recipe$var_info$role == "predictor"]
                        }, error = function(e) {
                                names(recipe$template)[names(recipe$template) != "value"]
                        })
                        
                        found <- intersect(required_preds, csv_cols)
                        missing <- setdiff(required_preds, csv_cols)
                        
                        div(
                                class = "imputed-summary mt-3",
                                div(class = "imputed-header",
                                    icon("table-columns"),
                                    sprintf(" %d / %d predictors found in CSV",
                                            length(found), length(required_preds))
                                ),
                                if (length(missing) > 0) {
                                        div(
                                                div(
                                                        class = "small mt-1 mb-1",
                                                        style = "color: var(--water-blue);",
                                                        sprintf("%d missing — will be imputed:", 
                                                                length(missing))
                                                ),
                                                div(
                                                        lapply(missing, function(v) {
                                                                label <- if (v %in% names(PREDICTOR_SPEC)) {
                                                                        PREDICTOR_SPEC[[v]]$label
                                                                } else v
                                                                span(class = "imputed-var", label)
                                                        })
                                                )
                                        )
                                } else {
                                        div(
                                                class = "small mt-1",
                                                style = "color: var(--accent);",
                                                icon("check-circle"),
                                                " All required predictors present"
                                        )
                                }
                        )
                }, error = function(e) NULL)
        })
        
        
        # Batch prediction with imputation ----------------------------------------
        
        observeEvent(input$batch_predict, {
                req(input$batch_file)
                req(input$batch_metric)
                req(input$batch_taxon_group)
                
                model_filename <- build_remote_filename(
                        input$batch_taxon_group, input$batch_metric, type = "model"
                )
                model_file <- fetch_remote_file(model_filename)
                
                if (is.null(model_file)) {
                        showNotification("Could not fetch model from Zenodo", type = "error")
                        return()
                }
                
                withProgress(message = "Batch prediction", value = 0, {
                        
                        tryCatch({
                                incProgress(0.05, detail = "Reading CSV...")
                                data <- fread(input$batch_file$datapath)
                                data_df <- as.data.frame(data)
                                
                                incProgress(0.10, detail = "Loading model...")
                                bundled <- readRDS(model_file)
                                model <- unbundle(bundled)
                                ranger_model <- extract_fit_engine(model)
                                recipe <- extract_recipe(model)
                                
                                # Identify required predictors
                                required_preds <- tryCatch({
                                        recipe$var_info$variable[recipe$var_info$role == "predictor"]
                                }, error = function(e) {
                                        names(recipe$template)[names(recipe$template) != "value"]
                                })
                                
                                csv_cols <- names(data_df)
                                found_cols <- intersect(required_preds, csv_cols)
                                missing_cols <- setdiff(required_preds, csv_cols)
                                
                                # ── Load training data for imputation ──
                                incProgress(0.05, detail = "Loading training data...")
                                td <- NULL
                                td_filename <- build_remote_filename(input$batch_taxon_group, type = "training_data")
                                td_file <- fetch_remote_file(td_filename)
                                if (!is.null(td_file)) {
                                        td <- tryCatch({
                                                out <- readRDS(td_file)
                                                if (inherits(out, "data.table")) out <- as.data.frame(out)
                                                out
                                        }, error = function(e) NULL)
                                }
                                
                                impute_method <- input$batch_impute_method
                                k <- input$batch_knn_k %||% 50
                                
                                # ── Impute missing predictor columns ──
                                if (length(missing_cols) > 0) {
                                        n_imputed <- length(missing_cols)
                                        
                                        if (impute_method == "conditional" && !is.null(td)) {
                                                incProgress(0.05, 
                                                            detail = sprintf("Imputing %d columns via KNN...", 
                                                                             n_imputed))
                                                
                                                # Pre-compute training data stats once
                                                specified_vars_td <- intersect(found_cols, names(td))
                                                
                                                if (length(specified_vars_td) > 0) {
                                                        train_spec <- as.matrix(td[, specified_vars_td, drop = FALSE])
                                                        col_sds <- apply(train_spec, 2, sd, na.rm = TRUE)
                                                        col_means <- colMeans(train_spec, na.rm = TRUE)
                                                        col_sds[col_sds == 0] <- 1
                                                        train_scaled <- scale(train_spec, center = col_means, scale = col_sds)
                                                        
                                                        k_use <- min(k, max(10, floor(nrow(td) / 5)))
                                                        k_use <- min(k_use, nrow(td))
                                                        
                                                        unspec_td <- intersect(missing_cols, names(td))
                                                        
                                                        # Pre-allocate columns
                                                        for (mc in missing_cols) {
                                                                data_df[[mc]] <- NA_real_
                                                        }
                                                        
                                                        # Row-wise KNN with progress
                                                        n_rows <- nrow(data_df)
                                                        progress_chunk <- max(1, floor(n_rows / 20))
                                                        
                                                        for (i in seq_len(n_rows)) {
                                                                query_vals <- as.numeric(data_df[i, specified_vars_td])
                                                                query_scaled <- (query_vals - col_means) / col_sds
                                                                
                                                                diffs <- sweep(train_scaled, 2, query_scaled, "-")
                                                                dists <- sqrt(rowSums(diffs^2, na.rm = TRUE))
                                                                knn_idx <- order(dists)[1:k_use]
                                                                
                                                                for (mc in unspec_td) {
                                                                        data_df[i, mc] <- mean(td[knn_idx, mc], na.rm = TRUE)
                                                                }
                                                                
                                                                # Fallback for any missing_cols not in td
                                                                still_na <- setdiff(missing_cols, unspec_td)
                                                                for (mc in still_na) {
                                                                        data_df[i, mc] <- if (mc %in% names(PREDICTOR_SPEC)) {
                                                                                PREDICTOR_SPEC[[mc]]$value
                                                                        } else 0
                                                                }
                                                                
                                                                # Update progress bar periodically
                                                                if (i %% progress_chunk == 0 || i == n_rows) {
                                                                        incProgress(
                                                                                0.40 * (progress_chunk / n_rows),
                                                                                detail = sprintf(
                                                                                        "KNN imputation: row %s / %s",
                                                                                        format(i, big.mark = ","),
                                                                                        format(n_rows, big.mark = ",")
                                                                                )
                                                                        )
                                                                }
                                                        }
                                                } else {
                                                        # No overlapping columns — fall back to marginal
                                                        for (mc in missing_cols) {
                                                                data_df[[mc]] <- if (mc %in% names(td)) {
                                                                        median(td[[mc]], na.rm = TRUE)
                                                                } else if (mc %in% names(PREDICTOR_SPEC)) {
                                                                        PREDICTOR_SPEC[[mc]]$value
                                                                } else 0
                                                        }
                                                        incProgress(0.40, detail = "Columns imputed (marginal fallback)")
                                                }
                                                
                                        } else if (impute_method == "marginal" && !is.null(td)) {
                                                incProgress(0.05, detail = "Imputing via training medians...")
                                                for (mc in missing_cols) {
                                                        data_df[[mc]] <- if (mc %in% names(td)) {
                                                                median(td[[mc]], na.rm = TRUE)
                                                        } else if (mc %in% names(PREDICTOR_SPEC)) {
                                                                PREDICTOR_SPEC[[mc]]$value
                                                        } else 0
                                                }
                                                incProgress(0.40, detail = "Columns imputed")
                                                
                                        } else {
                                                incProgress(0.05, detail = "Applying defaults...")
                                                for (mc in missing_cols) {
                                                        data_df[[mc]] <- if (mc %in% names(PREDICTOR_SPEC)) {
                                                                PREDICTOR_SPEC[[mc]]$value
                                                        } else 0
                                                }
                                                incProgress(0.40, detail = "Defaults applied")
                                        }
                                        
                                        # ── Enforce compositional constraints row-wise ──
                                        incProgress(0.05, detail = "Enforcing compositional constraints...")
                                        for (grp in COMPOSITIONAL_GROUPS) {
                                                spec_in_grp <- intersect(grp, found_cols)
                                                unspec_in_grp <- intersect(grp, missing_cols)
                                                
                                                if (length(unspec_in_grp) > 0 && length(spec_in_grp) > 0) {
                                                        for (i in seq_len(nrow(data_df))) {
                                                                sum_specified <- sum(data_df[i, spec_in_grp], na.rm = TRUE)
                                                                remaining <- max(0, 1 - sum_specified)
                                                                
                                                                raw_vals <- as.numeric(data_df[i, unspec_in_grp])
                                                                sum_raw <- sum(raw_vals, na.rm = TRUE)
                                                                
                                                                if (sum_raw > 0) {
                                                                        scaled <- (raw_vals / sum_raw) * remaining
                                                                } else {
                                                                        scaled <- rep(remaining / length(unspec_in_grp),
                                                                                      length(unspec_in_grp))
                                                                }
                                                                data_df[i, unspec_in_grp] <- scaled
                                                        }
                                                }
                                        }
                                        
                                        # ── Enforce ordered spatial constraints ──
                                        for (oc in ORDERED_CONSTRAINTS) {
                                                chain    <- oc$chain
                                                interior <- oc$interior
                                                endpoints <- oc$endpoints
                                                
                                                chain_in_data <- intersect(chain, names(data_df))
                                                if (length(chain_in_data) < 2) next
                                                
                                                lo_var <- endpoints[1]
                                                hi_var <- endpoints[2]
                                                imputed_interior  <- intersect(interior, missing_cols)
                                                specified_interior <- intersect(interior, found_cols)
                                                
                                                for (i in seq_len(nrow(data_df))) {
                                                        lo_val <- if (lo_var %in% names(data_df)) data_df[i, lo_var] else NA
                                                        hi_val <- if (hi_var %in% names(data_df)) data_df[i, hi_var] else NA
                                                        
                                                        # Midpoint or match logic
                                                        if (length(imputed_interior) == length(interior) &&
                                                            !is.na(lo_val) && !is.na(hi_val)) {
                                                                midpoint <- (lo_val + hi_val) / 2
                                                                for (v in imputed_interior) data_df[i, v] <- midpoint
                                                                
                                                        } else if (length(imputed_interior) > 0 && length(specified_interior) > 0) {
                                                                ref_val <- data_df[i, specified_interior[1]]
                                                                for (v in imputed_interior) data_df[i, v] <- ref_val
                                                        }
                                                        
                                                        # Safety-net clamping
                                                        vals <- as.numeric(data_df[i, chain_in_data])
                                                        for (j in 2:length(chain_in_data)) {
                                                                if (vals[j] < vals[j - 1] && chain_in_data[j] %in% missing_cols)
                                                                        vals[j] <- vals[j - 1]
                                                        }
                                                        for (j in (length(chain_in_data) - 1):1) {
                                                                if (vals[j] > vals[j + 1] && chain_in_data[j] %in% missing_cols)
                                                                        vals[j] <- vals[j + 1]
                                                        }
                                                        data_df[i, chain_in_data] <- vals
                                                }
                                        }
                                        
                                        # ── Round integer columns ──
                                        for (mc in missing_cols) {
                                                if (isTRUE(PREDICTOR_SPEC[[mc]]$is_integer)) {
                                                        data_df[[mc]] <- round(data_df[[mc]])
                                                }
                                        }
                                        
                                } else {
                                        incProgress(0.50, detail = "All predictor columns present")
                                }
                                
                                # ── Add any non-predictor columns the recipe expects ──
                                # The recipe may require ID/metadata columns (e.g. scheme_id,
                                # taxon, metric, year) that are not predictors but were present
                                # during training. Add them as NA so bake() doesn't fail.
                                incProgress(0.05, detail = "Preparing data for model...")
                                all_recipe_vars <- recipe$var_info$variable
                                still_missing <- setdiff(all_recipe_vars, names(data_df))
                                if (length(still_missing) > 0) {
                                        cat("  Adding", length(still_missing), 
                                            "non-predictor recipe columns as NA:",
                                            paste(still_missing, collapse = ", "), "\n")
                                        for (mv in still_missing) {
                                                data_df[[mv]] <- NA
                                        }
                                }
                                
                                # ── Bake and predict ──
                                incProgress(0.10, detail = "Running predictions...")
                                prepped_data <- bake(recipe, new_data = as.data.table(data_df))
                                
                                quantiles <- sort(as.numeric(input$batch_quantiles))
                                if (length(quantiles) == 0) {
                                        quantiles <- c(0.1, 0.5, 0.9)
                                }
                                qpred <- predict(
                                        ranger_model,
                                        data = prepped_data,
                                        type = "quantiles",
                                        quantiles = quantiles
                                )
                                
                                # ── Assemble results ──
                                incProgress(0.05, detail = "Assembling results...")
                                # Use the original uploaded data (not the imputed frame) as the
                                # left-hand side so users see their own columns
                                pred_df <- as.data.frame(qpred$predictions)
                                pred_df[] <- lapply(pred_df, round, digits = 2)
                                
                                results <- cbind(data, pred_df)
                                names(results)[(ncol(data)+1):ncol(results)] <- paste0("Q", quantiles * 100)
                                
                                batch_results(results)
                                
                                incProgress(0.10, detail = "Done!")
                                showNotification(
                                        paste("Predictions complete for", nrow(data), "rows"),
                                        type = "message"
                                )
                                
                        }, error = function(e) {
                                showNotification(paste("Error:", e$message), type = "error")
                        })
                        
                }) # end withProgress
        })
        
        output$batch_results <- renderDT({
                req(batch_results())
                datatable(
                        batch_results(),
                        options = list(
                                pageLength = 10,
                                scrollX = TRUE,
                                dom = 'Bfrtip'
                        ),
                        class = 'cell-border stripe'
                )
        })
        
        output$download_results <- downloadHandler(
                filename = function() {
                        paste0("qrf_predictions_", input$batch_taxon_group, "_", 
                               input$batch_metric, "_", Sys.Date(), ".csv")
                },
                content = function(file) {
                        req(batch_results())
                        fwrite(batch_results(), file)
                }
        )
        
        # Parameter reference table
        output$param_reference <- renderDT({
                ref_data <- data.frame(
                        Parameter = names(PREDICTOR_SPEC),
                        Label = sapply(PREDICTOR_SPEC, function(x) x$label),
                        Group = sapply(PREDICTOR_SPEC, function(x) x$group),
                        Type = sapply(PREDICTOR_SPEC, function(x) x$type),
                        Default = sapply(PREDICTOR_SPEC, function(x) as.character(x$value))
                )
                
                datatable(
                        ref_data,
                        options = list(
                                pageLength = 25,
                                dom = 't',
                                ordering = FALSE
                        ),
                        rownames = FALSE,
                        class = 'compact'
                )
        })
}



# Run App


shinyApp(ui = ui, server = server)