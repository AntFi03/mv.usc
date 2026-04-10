load_all()
library(tidyverse)


setwd(
  "/home/anton/Documentos/MTE Programación/C2/Análise multivariante/Traballos/Traballo 1"
)
dat_source <- read.table("ingenieros_pilotos.txt", header = TRUE) |>
  rownames_to_column(var = "id")
head(dat_source)


# --- Contraste De Igualdade Das Matrices De Covarianza -- #
test_compare_covariance(
  dat_list = list(as.matrix(dat_source[, 2:7]), as.matrix(dat_source[, 8:13]))
)
test_compare_covariance()
