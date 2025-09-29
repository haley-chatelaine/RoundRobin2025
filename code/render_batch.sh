#!/bin/bash

# Loop through all Rmd files in the current directory
for rmd_file in *.Rmd; do
  echo "Rendering $rmd_file..."
  Rscript -e "rmarkdown::render('$rmd_file')"
  if [ $? -eq 0 ]; then
    echo "$rmd_file rendered successfully."
  else
    echo "Error rendering $rmd_file."
  fi
done