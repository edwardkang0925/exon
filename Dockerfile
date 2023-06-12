FROM r-base:4.2.2

ENV RENV_VERSION 0.17.3
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    software-properties-common \
    dirmngr \
    wget \
    build-essential \
    libssl-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libtiff5-dev \
    cmake

# Clean up
RUN apt-get autoremove -y

WORKDIR /project
COPY renv.lock renv.lock

RUN mkdir -p renv
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
#COPY renv/settings.dcf renv/settings.dcf

# note: update this path as necessary based no the r-base r version
# and what you make your WORKDIR
# x86_64-pc-linux-gnu for linux
ENV R_LIBS /project/renv/library/R-4.2/x86_64-pc-linux-gnu

RUN R -e "renv::restore()"
