FROM rocker/r-ver:4.4.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
    libtiff5-dev libjpeg-dev libhdf5-dev make pandoc git && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /work
COPY . /work
RUN Rscript scripts/install_packages.R
CMD ["Rscript", "scripts/run_pipeline.R"]
