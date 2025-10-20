#' dist.pois.nbinom.hamim: GOF for Poisson & Negative Binomial with auto pooling
#'
#' **EN:** Tools to preview and test χ² GOF for Poisson/Negative-Binomial on discrete counts,
#' with two-sided automatic pooling and an optional single tail bin ("t+").
#'
#' **ID:** Alat untuk pratinjau dan uji χ² GOF untuk Poisson/Binomial Negatif pada data cacah,
#' dengan pooling dua-arah otomatis dan opsi satu kelas ekor ("t+").
#'
#' @section Quick Start (EN):
#' \preformatted{
#' library(dist.pois.nbinom.hamim)
#' txt <- "X Fx\n0 114\n1 25\n2 15\n3 10\n4 6\n5 5\n6 2\n7 1\n8 1\n9 0\n10 1\n"
#' dat <- read.table(text = txt, header = TRUE)
#' prev <- distribusi.hamim(dat, "poisson", pool="none", tail=TRUE, plot=FALSE)
#' prev$tabel
#' }
#'
#' @section Mulai Cepat (ID):
#' \preformatted{
#' library(dist.pois.nbinom.hamim)
#' tempel <- "X Fx\n0 114\n1 25\n2 15\n3 10\n4 6\n5 5\n6 2\n7 1\n8 1\n9 0\n10 1\n"
#' dataku <- read.table(text = tempel, header = TRUE)
#' uji_nb <- distribusi.hamim(dataku, "negbinom", method="log.iterasi", ex_min=5, pool="auto", tail=TRUE)
#' uji_nb$p_val
#' }
#'
#' @keywords package
"_PACKAGE"
