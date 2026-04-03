parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) {
      stop("Arguments must use --key=value format.")
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- pieces[1]
    value <- if (length(pieces) > 1) paste(pieces[-1], collapse = "=") else "TRUE"
    out[[key]] <- value
  }
  out
}

as_bool <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  tolower(x) %in% c("true", "1", "t", "yes", "y")
}

make_beta_true <- function(p) {
  signal <- c(2.0, -1.8, 1.5, -1.2, 1.0, 0.8, -0.6, 0.5)
  beta_true <- numeric(p)
  beta_true[seq_len(min(length(signal), p))] <- signal[seq_len(min(length(signal), p))]
  beta_true
}

simulate_data <- function(n, p, sigma_eps, data_seed) {
  set.seed(data_seed)
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta_true <- make_beta_true(p)
  y <- drop(x %*% beta_true + rnorm(n, sd = sigma_eps))
  list(X = x, Y = y, beta_true = beta_true)
}

seed_for_prior <- function(prior) {
  switch(prior,
    horseshoe = 1101L,
    laplace = 1102L,
    ridge = 1103L,
    sharkfin = 1104L,
    nonlocal = 1105L,
    inverselaplace = 1106L,
    stop("Unknown prior: ", prior)
  )
}

extra_prior_args <- function(prior, p) {
  if (prior == "sharkfin") {
    return(list(prob_vec = rep(0.25, p)))
  }
  if (prior == "nonlocal") {
    return(list(prior_mean = rep(1.5, p)))
  }
  if (prior == "inverselaplace") {
    return(list(lambda = 1.0))
  }
  list()
}

run_one_prior <- function(prior, lib, n, p, N, burnin, sigma_eps, standardize, singular, data_seed) {
  library(bayeslm, lib.loc = lib)
  dat <- simulate_data(n = n, p = p, sigma_eps = sigma_eps, data_seed = data_seed)

  args <- c(
    list(
      Y = dat$Y,
      X = dat$X,
      prior = prior,
      penalize = rep(1, p),
      block_vec = rep(1, p),
      sigma = NULL,
      N = N,
      burnin = burnin,
      thinning = 1L,
      vglobal = 1,
      sampling_vglobal = TRUE,
      verb = FALSE,
      icept = FALSE,
      standardize = standardize,
      singular = singular,
      scale_sigma_prior = TRUE,
      cc = rep(1, p)
    ),
    extra_prior_args(prior, p)
  )

  set.seed(seed_for_prior(prior))
  elapsed <- unname(system.time(fit <- do.call(bayeslm, args))[["elapsed"]])
  beta_mean <- colMeans(fit$beta)

  out <- list(
    prior = prior,
    elapsed = elapsed,
    beta_mean = beta_mean,
    beta_err_max = max(abs(beta_mean - dat$beta_true)),
    beta_rmse = sqrt(mean((beta_mean - dat$beta_true)^2)),
    sigma_mean = mean(fit$sigma)
  )

  if (!is.null(fit$vglobal)) {
    out$vglobal_mean <- mean(fit$vglobal)
  }
  if (!is.null(fit$tau)) {
    out$tau_mean <- mean(fit$tau)
  }
  if (!is.null(fit$lambda)) {
    out$lambda_mean <- mean(fit$lambda)
  }

  out
}

run_child <- function(opts) {
  required <- c("lib", "outfile", "n", "p", "N", "burnin", "sigma_eps", "standardize", "priors", "data_seed")
  missing <- required[!required %in% names(opts)]
  if (length(missing) > 0) {
    stop("Missing required args: ", paste(missing, collapse = ", "))
  }

  if (!is.null(opts$threads)) {
    suppressPackageStartupMessages(library(RcppParallel))
    RcppParallel::setThreadOptions(numThreads = as.integer(opts$threads))
  }

  priors <- strsplit(opts$priors, ",", fixed = TRUE)[[1]]
  priors <- priors[nzchar(priors)]
  out <- list(
    meta = list(
      lib = opts$lib,
      n = as.integer(opts$n),
      p = as.integer(opts$p),
      N = as.integer(opts$N),
      burnin = as.integer(opts$burnin),
      sigma_eps = as.numeric(opts$sigma_eps),
      standardize = as_bool(opts$standardize),
      singular = as_bool(opts$singular),
      data_seed = as.integer(opts$data_seed),
      priors = priors,
      beta_true = make_beta_true(as.integer(opts$p))
    )
  )

  for (prior in priors) {
    out[[prior]] <- run_one_prior(
      prior = prior,
      lib = opts$lib,
      n = as.integer(opts$n),
      p = as.integer(opts$p),
      N = as.integer(opts$N),
      burnin = as.integer(opts$burnin),
      sigma_eps = as.numeric(opts$sigma_eps),
      standardize = as_bool(opts$standardize),
      singular = as_bool(opts$singular),
      data_seed = as.integer(opts$data_seed)
    )
  }

  saveRDS(out, opts$outfile)
}

metric_gap <- function(candidate, benchmark, metric) {
  if (is.null(candidate[[metric]]) || is.null(benchmark[[metric]])) {
    return(NA_real_)
  }
  candidate[[metric]] - benchmark[[metric]]
}

run_driver <- function(opts) {
  required <- c("benchmark_lib", "candidate_lib")
  missing <- required[!required %in% names(opts)]
  if (length(missing) > 0) {
    stop("Missing required args: ", paste(missing, collapse = ", "))
  }

  priors <- if (!is.null(opts$priors)) opts$priors else "horseshoe,laplace,ridge,sharkfin,nonlocal,inverselaplace"
  n <- if (!is.null(opts$n)) as.integer(opts$n) else 600L
  p <- if (!is.null(opts$p)) as.integer(opts$p) else 120L
  N <- if (!is.null(opts$N)) as.integer(opts$N) else 150L
  burnin <- if (!is.null(opts$burnin)) as.integer(opts$burnin) else 50L
  sigma_eps <- if (!is.null(opts$sigma_eps)) as.numeric(opts$sigma_eps) else 0.1
  standardize <- if (!is.null(opts$standardize)) opts$standardize else "TRUE"
  singular <- if (!is.null(opts$singular)) opts$singular else "FALSE"
  data_seed <- if (!is.null(opts$data_seed)) as.integer(opts$data_seed) else 20260430L
  threads <- if (!is.null(opts$threads)) as.integer(opts$threads) else 1L
  outfile <- if (!is.null(opts$outfile)) opts$outfile else tempfile(fileext = ".rds")

  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) != 1) {
    stop("Could not determine script path from commandArgs().")
  }
  script <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
  bench_file <- tempfile(fileext = ".rds")
  cand_file <- tempfile(fileext = ".rds")

  # Load the two package builds in separate R processes so each run gets a
  # clean namespace. Reusing one process can silently compare the same loaded
  # DLL twice, which is especially misleading for exact MCMC regression tests.
  common_args <- c(
    paste0("--mode=child"),
    paste0("--n=", n),
    paste0("--p=", p),
    paste0("--N=", N),
    paste0("--burnin=", burnin),
    paste0("--sigma_eps=", sigma_eps),
    paste0("--standardize=", standardize),
    paste0("--singular=", singular),
    paste0("--priors=", priors),
    paste0("--data_seed=", data_seed)
  )
  if (!is.na(threads)) {
    common_args <- c(common_args, paste0("--threads=", threads))
  }

  bench_status <- system2(
    command = file.path(R.home("bin"), "Rscript"),
    args = c(script, common_args, paste0("--lib=", opts$benchmark_lib), paste0("--outfile=", bench_file))
  )
  if (bench_status != 0) {
    stop("Benchmark run failed.")
  }

  cand_status <- system2(
    command = file.path(R.home("bin"), "Rscript"),
    args = c(script, common_args, paste0("--lib=", opts$candidate_lib), paste0("--outfile=", cand_file))
  )
  if (cand_status != 0) {
    stop("Candidate run failed.")
  }

  benchmark <- readRDS(bench_file)
  candidate <- readRDS(cand_file)
  prior_names <- benchmark$meta$priors

  summary_rows <- lapply(prior_names, function(prior) {
    bench <- benchmark[[prior]]
    cand <- candidate[[prior]]
    data.frame(
      prior = prior,
      benchmark_elapsed = bench$elapsed,
      candidate_elapsed = cand$elapsed,
      speedup = bench$elapsed / cand$elapsed,
      benchmark_beta_err_max = bench$beta_err_max,
      candidate_beta_err_max = cand$beta_err_max,
      beta_err_gap = cand$beta_err_max - bench$beta_err_max,
      benchmark_beta_rmse = bench$beta_rmse,
      candidate_beta_rmse = cand$beta_rmse,
      beta_rmse_gap = cand$beta_rmse - bench$beta_rmse,
      beta_mean_maxdiff = max(abs(cand$beta_mean - bench$beta_mean)),
      beta_mean_rmsdiff = sqrt(mean((cand$beta_mean - bench$beta_mean)^2)),
      sigma_gap = cand$sigma_mean - bench$sigma_mean,
      vglobal_gap = metric_gap(cand, bench, "vglobal_mean"),
      tau_gap = metric_gap(cand, bench, "tau_mean"),
      lambda_gap = metric_gap(cand, bench, "lambda_mean")
    )
  })

  summary_table <- do.call(rbind, summary_rows)
  result <- list(
    meta = benchmark$meta,
    benchmark = benchmark,
    candidate = candidate,
    summary = summary_table
  )

  saveRDS(result, outfile)
  print(summary_table, row.names = FALSE, digits = 6)
  cat("\nSaved full result to", outfile, "\n")
}

"%||%" <- function(x, y) {
  if (is.null(x)) y else x
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
mode <- if (!is.null(opts$mode)) opts$mode else "driver"

if (mode == "child") {
  run_child(opts)
} else if (mode == "driver") {
  run_driver(opts)
} else {
  stop("Unknown mode: ", mode)
}
