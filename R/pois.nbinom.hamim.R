#' Info singkat indeks dispersi & GOF (EN/ID)
#' @rdname info.pois.nbinom.hamim
#' @export

info.distribusi.hamim <- function() {
  cat("\n===== Penjelasan Statistik Distribusi =====\n")
  cat("1. Index of Dispersion (ID):\n")
  cat("   - Rasio antara varians dan mean.\n")
  cat("   - ID > 1 \u2192 overdispersed (penyebaran lebih besar dari acak).\n")
  cat("   - ID \u2248 1 \u2192 pola acak (Poisson).\n")
  cat("\n2. Green's Index (GI):\n")
  cat("   - Ukuran clumping/spasial agregasi: GI = (ID - 1)/(N - 1).\n")
  cat("   - GI > 0 menunjukkan adanya pengelompokan (clumping).\n")
  cat("\n3. d_statistic:\n")
  cat("   - Transformasi dari chi-square: d = sqrt(2*X^2) - sqrt(2N - 3).\n")
  cat("   - |d| > 1.96 \u2192 hasil signifikan pada \u03B1 = 0.05.\n")
  cat("\n4. Chi-square goodness-of-fit:\n")
  cat("   - Digunakan untuk menguji kecocokan data dengan distribusi tertentu.\n")
  cat("   - Jika p-value > 0.05, data sesuai dengan distribusi.\n")
  cat("\n====================================================\n")
  invisible(NULL)
}

# Alias nama baru agar keduanya bisa dipakai
#' @rdname info.pois.nbinom.hamim
#' @export
info.pois.nbinom.hamim <- info.distribusi.hamim


# ===== Helper: pooling dua-arah, dengan opsi kunci ekor =====
.pool_gof_bins <- function(df, ex_min = 5, lock_tail = FALSE, lock_head = FALSE) {
  stopifnot(all(c("num_X", "Fx", "Px", "Ex") %in% names(df)))
  groups <- data.frame(
    lo = df$num_X, hi = df$num_X,
    Fx = df$Fx, Px = df$Px, Ex = df$Ex,
    stringsAsFactors = FALSE
  )
  while (any(groups$Ex < ex_min) && nrow(groups) > 1) {
    # pilih kandidat indeks dengan Ex minimum, hormati kunci head/tail
    idx <- which(groups$Ex < ex_min)
    if (lock_head) idx <- idx[idx != 1]
    if (lock_tail) idx <- idx[idx != nrow(groups)]
    # jika semua sel di bawah ambang adalah head/tail yang terkunci,
    # tetap izinkan gabung, tapi arah gabung dipaksa sesuai kunci
    if (length(idx) == 0) idx <- which.min(groups$Ex)
    i <- idx[which.min(groups$Ex[idx])]
    can_left  <- i > 1
    can_right <- i < nrow(groups)
    # jika i adalah head terkunci, hanya boleh ke kanan
    if (lock_head && i == 1) can_left <- FALSE
    # jika i adalah tail terkunci, hanya boleh ke kiri
    if (lock_tail && i == nrow(groups)) can_right <- FALSE

    if (!can_left && !can_right) break
    if (can_left && can_right) {
      ExL <- groups$Ex[i-1] + groups$Ex[i]
      ExR <- groups$Ex[i]   + groups$Ex[i+1]
      dir <- if (ExR >= ExL) "right" else "left"
    } else if (can_left) {
      dir <- "left"
    } else {
      dir <- "right"
    }

    if (dir == "left") {
      j <- i - 1
      new <- data.frame(
        lo = groups$lo[j],
        hi = groups$hi[i],
        Fx = groups$Fx[j] + groups$Fx[i],
        Px = groups$Px[j] + groups$Px[i],
        Ex = groups$Ex[j] + groups$Ex[i]
      )
      groups <- rbind(
        if (j > 1) groups[1:(j-1), ] else NULL,
        new,
        if (i < nrow(groups)) groups[(i+1):nrow(groups), ] else NULL
      )
    } else {
      j <- i + 1
      new <- data.frame(
        lo = groups$lo[i],
        hi = groups$hi[j],
        Fx = groups$Fx[i] + groups$Fx[j],
        Px = groups$Px[i] + groups$Px[j],
        Ex = groups$Ex[i] + groups$Ex[j]
      )
      groups <- rbind(
        if (i > 1) groups[1:(i-1), ] else NULL,
        new,
        if (j < nrow(groups)) groups[(j+1):nrow(groups), ] else NULL
      )
    }
    rownames(groups) <- NULL
  }
  lab <- ifelse(groups$lo == groups$hi,
                as.character(groups$lo),
                paste0(groups$lo, "\u2013", groups$hi))
  data.frame(X = lab, Fx = groups$Fx, Px = groups$Px, Ex = groups$Ex, stringsAsFactors = FALSE)
}

#' Goodness-of-fit for Poisson / NegBin with two-sided pooling (EN/ID)
#'
#' @param data data.frame dengan kolom X (cacah, integer) dan Fx (frekuensi).
#' @param distrib c("poisson","negbinom")
#' @param method c("MLE","k.langsung","manual.iterasi","log.iterasi")
#' @param ex_min numeric; ambang Ex
#' @param min_expected alias lama untuk ex_min
#' @param plot,tail,tail_from,pool lihat help
#'
#' @return list: estimasi & tabel hasil
#'
#' @examples
#' dat <- data.frame(X = 0:6, Fx = c(114, 25, 15, 10, 6, 5, 2))
#' prev <- distribusi.hamim(dat, distrib = "poisson", pool = "none", tail = TRUE, plot = FALSE)
#' sum(prev$tabel$Px)
#' gof  <- distribusi.hamim(dat, distrib = "poisson", ex_min = 3, pool = "auto", tail = TRUE, plot = FALSE)
#' gof$p_val; gof$df
#' nb   <- distribusi.hamim(dat, distrib = "negbinom", method = "log.iterasi",
#'                          ex_min = 5, pool = "auto", tail = TRUE, plot = FALSE)
#' nb$k; nb$mu; nb$p_val
#' @export
# ====== FUNGSI UTAMA (bagian yang berubah ditandai komentar 'BARU') ======
distribusi.hamim <- function(data,
                             distrib = c("poisson", "negbinom"),
                             method = c("MLE", "k.langsung", "manual.iterasi", "log.iterasi"),
                             ex_min = 5,
                             min_expected = NULL,
                             plot = TRUE,
                             tail = FALSE,
                             tail_from = NULL,
                             pool = c("auto","none")) {  # [BARU]
  distrib <- match.arg(distrib)
  method  <- match.arg(method)
  pool    <- match.arg(pool)       # [BARU]
  if (!is.null(min_expected)) ex_min <- min_expected
  stopifnot(all(c("X","Fx") %in% names(data)))

  # --- siapkan data
  data$X <- as.numeric(as.character(data$X))
  data <- data[order(data$X), ]
  rownames(data) <- NULL
  x  <- data$X
  Fx <- data$Fx

  N <- sum(Fx)
  n <- sum(x * Fx)
  mean_x <- n / N
  var_x  <- if (N > 1) sum((x - mean_x)^2 * Fx) / (N - 1) else 0

  ID <- if (mean_x > 0) var_x / mean_x else NA_real_
  GI <- if (!is.na(ID) && N > 1) (ID - 1) / (N - 1) else NA_real_
  chi_ID <- if (!is.na(ID) && N > 1) ID * (N - 1) else NA_real_
  p_chi_ID <- if (!is.na(chi_ID)) stats::pchisq(chi_ID, df = N - 1, lower.tail = FALSE) else NA_real_
  d_stat <- if (!is.na(chi_ID)) sqrt(2 * chi_ID) - sqrt(2 * (N - 1) - 1) else NA_real_

  k0 <- if (!is.na(var_x) && var_x > mean_x) (mean_x^2) / (var_x - mean_x) else NA_real_

  # --- Estimasi parameter & PMF
  method_used <- ""
  k <- NA_real_; mu <- NA_real_; pmf <- NULL

  if (distrib == "poisson") {
    mu <- mean_x
    pmf <- function(xx) stats::dpois(xx, lambda = mu)
    method_used <- "Poisson: lambda = mean"
  } else {
    if (method == "MLE") {
      if (!requireNamespace("MASS", quietly = TRUE)) stop("Paket 'MASS' diperlukan untuk MLE.")
      obs_data <- rep(x, Fx)
      fit <- try(suppressWarnings(MASS::fitdistr(obs_data, "Negative Binomial")), silent = TRUE)
      if (inherits(fit, "try-error")) {
        if (var_x <= mean_x) stop("MLE gagal dan var <= mean; data tidak cocok untuk NegBin.")
        k <- (mean_x^2) / (var_x - mean_x); mu <- mean_x
        method_used <- "NegBin: fallback k.langsung (MLE gagal)"
      } else {
        k  <- as.numeric(fit$estimate["size"])
        mu <- as.numeric(fit$estimate["mu"])
        method_used <- "NegBin: Maximum Likelihood Estimation (MLE)"
      }
    } else if (method == "k.langsung") {
      if (var_x <= mean_x) stop("Metode k.langsung tidak cocok karena var <= mean.")
      k <- (mean_x^2) / (var_x - mean_x); mu <- mean_x
      method_used <- "NegBin: k = (mean^2)/(var - mean)"
    } else if (method == "manual.iterasi") {
      if (!any(x == 0)) stop("manual.iterasi membutuhkan kelas X == 0.")
      lhs <- (N - Fx[x == 0]) / N
      ks <- seq(0.05, 20, by = 0.0005)
      rhs <- ks * (mean_x / (ks + mean_x))
      k <- ks[which.min(abs(rhs - lhs))]; mu <- mean_x
      method_used <- "NegBin: Manual Iterasi berdasar proporsi nol"
    } else if (method == "log.iterasi") {
      if (!any(x == 0)) stop("log.iterasi membutuhkan kelas X == 0.")
      F0 <- Fx[x == 0]; if (length(F0) == 0 || F0 == 0) stop("log.iterasi membutuhkan F0 > 0.")
      lhs <- log10(N / F0)
      f_k <- function(k) abs(k * log10(1 + mean_x / k) - lhs)
      k  <- stats::optimize(f_k, interval = c(0.01, 100))$minimum
      mu <- mean_x
      method_used <- "NegBin: Log Iterasi (berbasis proporsi nol)"
    }
    pmf <- function(xx) stats::dnbinom(xx, size = k, mu = mu)
  }

  # CDF sesuai distribusi (dipakai untuk hitung massa ekor)
  cdf <- function(xx) {
    if (distrib == "poisson") stats::ppois(xx, lambda = mu)
    else                      stats::pnbinom(xx, size = k, mu = mu)
  }

  # --- Bangun tabel dasar: mode PRATINJAU (tanpa pooling) atau AUTO
  #     Agar "seluruh kelas" tercakup saat preview, tail default = (max(X)+1)+
  build_df0 <- function() {
    if (!isTRUE(tail)) {
      # hanya kelas yang ada di data (massa peluang bisa < 1)
      Px <- pmf(x); Ex <- Px * N
      return(data.frame(num_X = x, Fx = Fx, Px = Px, Ex = Ex,
                        stringsAsFactors = FALSE))
    }
    # tail = TRUE:
    if (pool == "none" && is.null(tail_from)) {
      t_from <- max(x) + 1L    # default pratinjau: seluruh ekor > max(X) jadi satu kelas
    } else if (pool == "auto" && is.null(tail_from)) {
      # 1) Cari kelas pertama yang memiliki Ex_x < ex_min
      find_first_below <- function() {
        xk <- 0L
        repeat {
          if (N * pmf(xk) < ex_min) return(xk)
          xk <- xk + 1L
          if (xk > 1e5) return(xk)  # guard
        }
      }
      t <- find_first_below()  # contoh Poisson webworm: t = 5

      # 2) Geser ke kiri sampai massa ekor memenuhi ambang (Ex_tail >= ex_min)
      repeat {
        Ex_tail <- N * (1 - cdf(t - 1))
        if (Ex_tail >= ex_min || t <= 1L) break
        t <- t - 1L
      }
      t_from <- t  # contoh: 4 -> kelas ekor jadi "4+"
    } else {



      t_from <- as.integer(tail_from)
      if (t_from < 1) stop("tail_from harus >= 1.")
    }

    base_x  <- 0:(t_from - 1)
    Px_head <- pmf(base_x)
    Px_tail <- max(0, 1 - sum(Px_head))    # proteksi numerik
    Ex_head <- N * Px_head
    Ex_tail <- N * Px_tail

    Fx_head <- sapply(base_x, function(xx) { idx <- which(x == xx); if (length(idx)) Fx[idx] else 0 })
    Fx_tail <- sum(Fx[x >= t_from])

    df <- data.frame(
      num_X = c(base_x, t_from),
      Fx    = c(Fx_head, Fx_tail),
      Px    = c(Px_head, Px_tail),
      Ex    = c(Ex_head, Ex_tail),
      stringsAsFactors = FALSE
    )
    attr(df, "t_from") <- t_from
    df
  }

  df0 <- if (isTRUE(tail) || pool == "none") build_df0() else {
    # tidak pakai tail sama sekali
    Px <- pmf(x); Ex <- Px * N
    data.frame(num_X = x, Fx = Fx, Px = Px, Ex = Ex, stringsAsFactors = FALSE)
  }
  t_from <- attr(df0, "t_from", exact = TRUE)

  # --- Jika pool="none": langsung tampilkan tanpa pooling & tanpa GOF
  if (pool == "none") {
    hasil_final <- within(df0, {
      X <- as.character(num_X)
      if (!is.null(t_from)) {
        X[num_X == t_from] <- paste0(t_from, "+")
      }
    })
    hasil_final <- hasil_final[, c("X","Fx","Px","Ex")]
    chi_sq <- NA_real_; df_chi <- NA_integer_; p_val <- NA_real_

    cat("\n================= PRATINJAU (tanpa pooling) =================\n")
    cat("Distribusi yang diuji     :", distrib, "\n")
    cat("Metode estimasi           :", method_used, "\n")
    if (!is.na(k0)) cat("Dugaan awal k (k0)         :", round(k0,4), "\n")
    if (!is.na(k))  cat("Nilai final k              :", round(k,4), "\n")
    cat("N, n, mean, var           :", N, ",", n, ",", round(mean_x,4), ",", round(var_x,4), "\n")
    cat("Catatan                   : Ini hanya pratinjau Px & Ex per kelas, belum uji \\u03C7^2.\n")
    cat("==============================================================\n")
    print(data.frame(
      X = hasil_final$X,
      Fx = hasil_final$Fx,
      Px = round(hasil_final$Px, 5),
      Ex = round(hasil_final$Ex, 2),
      check.names = FALSE
    ))

    if (isTRUE(plot)) {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("Paket 'ggplot2' tidak tersedia; plot dilewati.")
      } else {
        # set urutan level SEKALI di luar aes()
        hasil_final$X <- factor(hasil_final$X, levels = hasil_final$X)

        p <- ggplot2::ggplot(hasil_final, ggplot2::aes(x = X)) +  # pakai 'X' saja
          ggplot2::geom_bar(ggplot2::aes(y = Fx, fill = "Observed"), stat = "identity") +
          ggplot2::geom_point(ggplot2::aes(y = Ex, color = "Expected"), size = 3) +
          ggplot2::scale_fill_manual(name = "", values = c("Observed" = "skyblue")) +
          ggplot2::scale_color_manual(name = "", values = c("Expected" = "red")) +
          ggplot2::labs(
            title = paste0("Distribusi ", distrib,
                           if (pool == "none") ": Observed vs Expected (tanpa pooling)"
                           else ": Observed vs Expected (setelah pooling)"),
            x = "Kelas X", y = "Frekuensi"
          ) +
          ggplot2::theme_minimal()
        print(p)


      }
    }

    return(invisible(list(
      distrib = distrib, method = method_used, k0 = k0, k = k, mu = mu,
      N = N, n = n, mean = mean_x, varians = var_x, ID = ID, GI = GI,
      chi_ID = chi_ID, p_ID = p_chi_ID, d_stat = d_stat,
      chi_sq = chi_sq, df = df_chi, p_val = p_val,
      tabel = hasil_final,
      t_from = t_from
    )))
  }

  # --- Jika pool="auto": lanjut pooling dua-arah (seperti biasa)
  hasil_final <- .pool_gof_bins(df0, ex_min = ex_min, lock_tail = isTRUE(tail), lock_head = FALSE)
  if (isTRUE(tail) && !is.null(t_from)) {
    hasil_final$X <- sub(paste0("^", t_from, "$"), paste0(t_from, "+"), hasil_final$X)
    hasil_final$X <- sub(paste0("^", t_from, "\u2013Inf$"), paste0(t_from, "+"), hasil_final$X)
  }

  hasil_final$Chi_Square <- (hasil_final$Fx - hasil_final$Ex)^2 / hasil_final$Ex
  chi_sq <- sum(hasil_final$Chi_Square)
  m_par <- if (distrib == "poisson") 1 else 2
  df_chi <- nrow(hasil_final) - 1 - m_par
  p_val <- if (df_chi > 0) stats::pchisq(chi_sq, df = df_chi, lower.tail = FALSE) else NA_real_

  cat("\n======================= HASIL ANALISIS =======================\n")
  cat("Distribusi yang diuji     :", distrib, "\n")
  cat("Metode estimasi           :", method_used, "\n")
  if (!is.na(k0)) cat("Dugaan awal k (k0)         :", round(k0, 4), "\n")
  if (!is.na(k))  cat("Nilai final k              :", round(k, 4), "\n")
  cat("N (unit)                  :", N, "\n")
  cat("n (total individu)        :", n, "\n")
  cat("Rata-rata (mean)          :", round(mean_x, 4), "\n")
  cat("Varians                   :", round(var_x, 4), "\n")
  cat("Index Dispersion (ID)     :", round(ID, 4), "\n")
  cat("Green's Index             :", round(GI, 4), "\n")
  cat("Chi-square untuk ID       :", round(chi_ID, 4), ", p-value:", ifelse(is.na(p_chi_ID),"NA",round(p_chi_ID,4)), "\n")
  cat("d statistic               :", round(d_stat, 4), ifelse(!is.na(d_stat) && abs(d_stat) > 1.96, " (signifikan)", ""), "\n")
  cat("\n--- Uji Kesesuaian Distribusi (setelah pooling) ---\n")
  cat("Ambang Ex minimum         :", ex_min, "\n")
  if (isTRUE(tail) && !is.null(t_from)) cat("Batas ekor                : ", t_from, "+\n", sep = "")
  cat("Jumlah kelas akhir        :", nrow(hasil_final), "\n")
  cat("Chi-Square total          :", round(chi_sq, 4), "\n")
  cat("Derajat bebas             :", df_chi, "\n")
  cat("P-value                   :", ifelse(is.na(p_val),"NA",round(p_val,4)), "\n")
  if (!is.na(p_val)) {
    if (p_val > 0.05) cat("Kesimpulan                : Gagal tolak H0 (data cocok dengan distribusi)\n")
    else              cat("Kesimpulan                : Tolak H0 (data tidak cocok dengan distribusi)\n")
  }
  cat("============================================================\n")
  print(data.frame(
    X = hasil_final$X, Fx = hasil_final$Fx,
    Ex = round(hasil_final$Ex, 2),
    Chi_Square = round(hasil_final$Chi_Square, 2),
    check.names = FALSE
  ))

  if (isTRUE(plot)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Paket 'ggplot2' tidak tersedia; plot dilewati.")
    } else {

      # set urutan level SEKALI di luar aes()
      hasil_final$X <- factor(hasil_final$X, levels = hasil_final$X)

      p <- ggplot2::ggplot(hasil_final, ggplot2::aes(x = X)) +  # pakai 'X' saja
        ggplot2::geom_bar(ggplot2::aes(y = Fx, fill = "Observed"), stat = "identity") +
        ggplot2::geom_point(ggplot2::aes(y = Ex, color = "Expected"), size = 3) +
        ggplot2::scale_fill_manual(name = "", values = c("Observed" = "skyblue")) +
        ggplot2::scale_color_manual(name = "", values = c("Expected" = "red")) +
        ggplot2::labs(
          title = paste0("Distribusi ", distrib,
                         if (pool == "none") ": Observed vs Expected (tanpa pooling)"
                         else ": Observed vs Expected (setelah pooling)"),
          x = "Kelas X", y = "Frekuensi"
        ) +
        ggplot2::theme_minimal()
      print(p)


    }
  }

  invisible(list(
    distrib = distrib, method = method_used, k0 = k0, k = k, mu = mu,
    N = N, n = n, mean = mean_x, varians = var_x, ID = ID, GI = GI,
    chi_ID = chi_ID, p_ID = p_chi_ID, d_stat = d_stat,
    chi_sq = chi_sq, df = df_chi, p_val = p_val, tabel = hasil_final,
    t_from = t_from
  ))
}
