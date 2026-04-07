load_all()

test_mean(
  n = 20,
  cov.matrix = matrix(c(49, 35, 35, 52), nrow = 2),
  x.bar = c(102, 98),
  mu0 = c(100, 100)
)

?test_mean

test_mean(
  dat = read.table(
    paste(
      "http://eamo.usc.es/pub/pateiro/files/AM_learnr/data/",
      "notas.txt",
      sep = ""
    ),
    header = TRUE
  ),
  # cov.matrix = matrix(c(49, 35, 35, 52), nrow = 2),
  # mu0 = c(100, 100)
)
