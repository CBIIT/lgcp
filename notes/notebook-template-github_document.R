#' ---
#' title: "LGCP Notebook template"
#' author: "Brian Capaldo"
#' date: "4/12/2019"
#' output:
#'   github_document:
#'     toc: true
#' ---

#' Document preamble
#' 
#' This is a spin document, full details are available [here](https://rmarkdown.rstudio.com/articles_report_from_r_script.html). This is the template version for use at the LGCP at the NIH and all are welcome to distribute and modify.
#' This document will run as a regular Rscript if called using `source()` or will render into a github flavored markdown document if called using `rmarkdown::render()`.
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