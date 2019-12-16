# EmbryoTimecourse2018
This repo contains data and analysis scripts for our Embryo Timecourse paper. 
There are two types of content here: scripts to get the processed data (`download`), and scripts that we have used to do our analyses (`analysis_scripts`).

To use the download script, please ensure you have `curl` installed, and run the download script in the `download` folder. You can choose which of our datasets you would like to download while running the script.

Please note that you can now use our Bioconductor package [MouseGastrulationData](https://bioconductor.org/packages/release/data/experiment/html/MouseGastrulationData.html) to access these data in a more efficient way, with processed and raw count matrices delivered directly into your *R* session.
Download times should also be much shorter if you use the package - and you can use e.g. `loomr` to save the data to an easily readable format if you want to analyse the data using python.
