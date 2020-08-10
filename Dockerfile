## Loads an R image, based on the rocker 'geospatial' image
FROM rocker/binder:3.6.0

## Copy repo into ${HOME}, make user own $HOME
USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}
USER ${NB_USER}

## Run install.R from root directory
#RUN if [ -f install.R ]; then R --quiet -f install.R; fi



RUN apt-get install -y --no-install-recommends \
           cowplot \
           littler \

RUN installGithub.R stan-dev/rstanarm
