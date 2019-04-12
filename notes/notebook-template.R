#' ---
#' title: "LGCP Notebook template"
#' author: "Brian Capaldo"
#' date: "4/12/2019"
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#'     number_sections: true
#'     code_folding: show
#' ---

#' Document preamble
#' 
#' This is a spin document, full details are available [here](https://rmarkdown.rstudio.com/articles_report_from_r_script.html). This is the template version for use at the LGCP at the NIH and all are welcome to distribute and modify.
#' This document will run as a regular Rscript if called using `source()` or will render into a full featured html document if called using `rmarkdown::render()`.
#' As of the time of this writing, outputting to github does not support `toc_float`, `number_sections`, or `code_folding` options.
#' To render for github, change `output: ...` to `output: github_document`, or look at `notebook-template-github_document` for the `github_document` version of the template.
#' For additional details, see [here](https://rmarkdown.rstudio.com/github_document_format.html)
#'
#'# Introduction
#'
#'Introduce the project and any pertinent background.
#'
#'# Methods
#'
#'Describe the underlying approach to analyzing the data
#'
#'# Results
#'
#'Analysis code and plots goes here.

#+ r example R block

library(edgeR)
library(tidyverse)
library(annotables)

#'# references