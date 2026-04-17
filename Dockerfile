FROM rocker/shiny:4.4.0

# System dependencies for your R packages
RUN echo 'Acquire::ForceIPv4 "true";' > /etc/apt/apt.conf.d/99force-ipv4 && \
    apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# R packages
RUN R -e "install.packages(c( \
    'shiny', 'bslib', 'workflows', 'recipes', 'parsnip', 'bundle', \
    'ranger', 'data.table', 'DT', 'ggplot2', 'curl' \
    ), repos='https://cloud.r-project.org/')"

RUN echo 'run_as shiny; \n\
server { \n\
  listen 8080 0.0.0.0; \n\
  location / { \n\
    app_dir /srv/shiny-server; \n\
    log_dir /var/log/shiny-server; \n\
  } \n\
}' > /etc/shiny-server/shiny-server.conf

COPY app.R /srv/shiny-server/app.R



# Copy app
COPY app.R /srv/shiny-server/app.R

# Rahti runs containers with an arbitrary non-root UID, so directories
# the app writes to must be group-writable by root (GID 0). This applies
# to the cache director
#RUN mkdir -p /srv/shiny-server/cache && \
#    chgrp -R 0 /srv/shiny-server && \
#    chmod -R g=u /srv/shiny-server
RUN mkdir -p /srv/shiny-server/cache /srv/shiny-server/app_cache && \
    chgrp -R 0 /srv/shiny-server && \
    chmod -R g=u /srv/shiny-server && \
    chmod -R 777 /srv/shiny-server/cache /srv/shiny-server/app_cache


# Shiny needs to listen on all interfaces, not just localhost
ENV SHINY_HOST=0.0.0.0
ENV SHINY_PORT=8080

EXPOSE 8080

CMD ["/usr/bin/shiny-server"]
