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

confirm "Download Atlas data? [y/N]" && lftp -c "open -u jmlabftp,HOBICAmeer6 ftp2.cruk.cam.ac.uk; get embryo_reviewers/atlas_data.tar.gz"
confirm "Download WT-chimera data? [y/N]" && lftp -c "open -u jmlabftp,HOBICAmeer6 ftp2.cruk.cam.ac.uk; get embryo_reviewers/chimera_wt_data.tar.gz"
confirm "Download Tal1-chimera data? [y/N]" && lftp -c "open -u jmlabftp,HOBICAmeer6 ftp2.cruk.cam.ac.uk; get embryo_reviewers/chimera_tal1_data.tar.gz"
confirm "Download singularity images? [y/N]" && lftp -c "open -u jmlabftp,HOBICAmeer6 ftp2.cruk.cam.ac.uk; get embryo_reviewers/singularity.tar.gz"