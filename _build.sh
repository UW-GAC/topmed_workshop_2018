#!/usr/bin/env bash

#Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"

export R_LIBS=/projects/topmed/working_code/analysis_pipeline_genesis2/R_library:/projects/resources/gactools/R_packages/library

echo "bookdown::render_book('index.Rmd', 'bookdown::gitbook')" > tmp.R
R -q --vanilla < tmp.R
rm tmp.R
