# PULSE benchmark app for CSC Rahti (OpenShift)
#
# Two build targets:
#   app     (default) the app only; models are read from a volume mounted at MODEL_DIR
#   bundle  the app plus the models in models/, one image layer per taxon
#
# Rahti's internal registry accepts images of at most 5 GiB. The reduced models
# (about 7.3 GB) push the bundle past that, so the bundle goes to an external
# registry (ghcr.io) and the app target to Rahti. See build_push.sh.
#
# Build context (this folder):
#   Dockerfile, .dockerignore, app.R
#   models/   <taxon>__<metric>.rds and <taxon>__<metric>_intervals.rds (bundle only)

ARG R_VERSION=4.5.1
FROM rocker/r-ver:${R_VERSION} AS app

# Binary packages from the Posit Package Manager snapshot that rocker pins for this
# R version, so package versions are reproducible. RANGER_VERSION can pin ranger to
# the version the models were fitted with (versions$ranger in <stem>_run_info.rds).
ARG RANGER_VERSION=""
RUN install2.r --error --skipinstalled --ncpus -1 \
        shiny bslib ranger recipes data.table DT ggplot2 remotes \
 && if [ -n "${RANGER_VERSION}" ]; then \
        Rscript -e "remotes::install_version('ranger', version = '${RANGER_VERSION}', upgrade = 'never')"; \
    fi \
 && rm -rf /tmp/downloaded_packages /tmp/Rtmp*

WORKDIR /srv/app
COPY app.R /srv/app/app.R

# OpenShift runs containers as an arbitrary UID in group 0: give group 0 the owner's
# permissions on everything the app reads.
RUN mkdir -p /srv/app/models /srv/models \
 && chgrp -R 0 /srv/app /srv/models \
 && chmod -R g=u /srv/app /srv/models

# Defaults; the Deployment can override them (oc set env).
ENV MODEL_DIR=/srv/models \
    QRF_CACHE=1 \
    QRF_THREADS=1 \
    PORT=3838 \
    HOME=/tmp

EXPOSE 3838
USER 1001

# Shiny directly, without Shiny Server (its run_as fails under arbitrary UIDs)
CMD ["Rscript", "-e", "shiny::runApp('/srv/app', host = '0.0.0.0', port = as.integer(Sys.getenv('PORT', '3838')), launch.browser = FALSE)"]


FROM app AS bundle

# One layer per taxon: smaller layers upload and pull more reliably.
COPY --chown=1001:0 models/diatoms__*       /srv/app/models/
COPY --chown=1001:0 models/fish__*          /srv/app/models/
COPY --chown=1001:0 models/invertebrates__* /srv/app/models/
COPY --chown=1001:0 models/macrophytes__*   /srv/app/models/

ENV MODEL_DIR=/srv/app/models
