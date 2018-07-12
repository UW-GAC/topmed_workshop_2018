#!/usr/bin/env bash

#Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"

export R_LIBS=/projects/topmed/working_code/analysis_pipeline_2.0.1/R_library

echo "bookdown::render_book('index.Rmd', 'bookdown::gitbook')" > tmp.R
R-3.5.1 -q --vanilla < tmp.R
rm tmp.R
