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

#below each line is the arrayexpress link
#download seems slower, but probably the server is more reliable.
confirm "Download Atlas data? [y/N]" && curl https://content.cruk.cam.ac.uk/jmlab/atlas_data.tar.gz > atlas_data.tar.gz
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6967/E-MTAB-6967.processed.1.zip
confirm "Download WT-chimera data? [y/N]" && curl https://content.cruk.cam.ac.uk/jmlab/chimera_wt_data.tar.gz > chimera_wt_data.tar.gz
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7324/E-MTAB-7324.processed.1.zip
confirm "Download Tal1-chimera data? [y/N]" && curl https://content.cruk.cam.ac.uk/jmlab/chimera_tal1_data.tar.gz > chimera_tal1_data.tar.gz
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7325/E-MTAB-7325.processed.1.zip
confirm "Download singularity images? [y/N]" && curl https://content.cruk.cam.ac.uk/jmlab/singularity.tar.gz > singularity.tar.gz
# no AE equivalent