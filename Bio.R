source("MMM_code.R")
library(data.table)
library(R.utils)
transpose = Rfast::transpose
dat = fread("../data/GSM5238385_ME11_50um.fragments.tsv.gz")
spac = fread("../data/GSM5238385_ME11_50um_spatial.tar.gz")
