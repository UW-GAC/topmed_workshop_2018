#!/usr/bin/env bash

#Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"

echo "bookdown::render_book('index.Rmd', 'bookdown::gitbook')" > tmp.R
R -q --vanilla < tmp.R
rm tmp.R
