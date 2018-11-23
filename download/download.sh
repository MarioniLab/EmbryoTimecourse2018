#!/bin/bash

confirm() {
    # call with a prompt string or use a default
    read -r -p "${1:-Are you sure? [y/N]} " response
    case "$response" in
        [yY][eE][sS]|[yY]) 
            true
            ;;
        *)
            false
            ;;
    esac
}

confirm "Download Atlas data? [y/N]" && curl https://content.cruk.cam.ac.uk/jmlab/atlas_data.tar.gz > atlas_data.tar.gz
confirm "Download WT-chimera data? [y/N]" && curl https://content.cruk.cam.ac.uk/jmlab/chimera_wt_data.tar.gz > chimera_wt_data.tar.gz
confirm "Download Tal1-chimera data? [y/N]" && lftp -c curl https://content.cruk.cam.ac.uk/jmlab/chimera_tal1_data.tar.gz > chimera_tal1_data.tar.gz
confirm "Download singularity images? [y/N]" && lftp -c curl https://content.cruk.cam.ac.uk/jmlab/singularity.tar.gz > singularity.tar.gz