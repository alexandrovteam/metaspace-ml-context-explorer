FROM rocker/shiny-verse:4.1.0

RUN apt-get update && apt-get install -y libcurl4-openssl-dev \
    libfftw3-dev libjpeg-dev libtiff-dev \
    libxml2-dev libssl-dev \
    libxml2-dev && rm -rf /var/lib/apt/lists/*

ADD ./ShinyApp/ /srv/shiny-server/

COPY ./packages.R packages.R

RUN Rscript packages.R

