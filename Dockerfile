FROM --platform=linux/amd64 rocker/shiny:latest

# Install system libraries required by R packages + Python for generate-ont-report.py
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libuv1-dev \
    git \
    python3 \
    && rm -rf /var/lib/apt/lists/*

# Copy the Python report-generation script and make it executable
COPY generate-ont-report.py /usr/local/bin/generate-ont-report.py
RUN chmod +x /usr/local/bin/generate-ont-report.py

WORKDIR /app

# Copy renv files first for Docker layer caching
COPY renv.lock .Rprofile ./
COPY renv/activate.R renv/

# Restore R environment using renv
RUN R -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); install.packages('renv'); renv::restore(rebuild = TRUE)"

# Copy remaining app source files
COPY . .

EXPOSE 3838

ENV SCRIPT_PATH=/usr/local/bin/generate-ont-report.py

CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]
