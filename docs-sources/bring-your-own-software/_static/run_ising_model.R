#!/usr/bin/env Rscript

# run_ising_model.R
# 2D Ising model via Metropolis Monte Carlo (base R only)
# Usage:
#   Rscript run_ising_model.R --N 100 --temps "1.5,2.0,2.269,2.5" --nsweeps 10000 --burnin 1000 --thin 10 --output ising_results.csv
#
# Arguments:
#  --N       lattice size (N x N)
#  --temps   comma-separated list of temperatures
#  --nsweeps total number of sweeps per temperature
#  --burnin  number of initial sweeps to discard
#  --thin    record every 'thin' sweeps after burnin
#  --output  output CSV file

# Parse command line args
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default=NULL) {
  key <- paste0("--", name)
  pos <- match(key, args)
  if (!is.na(pos) && length(args) >= pos+1) return(args[pos+1])
  default
}
N       <- as.integer(get_arg("N", 100))
temps   <- as.numeric(strsplit(get_arg("temps", "2.269"), ",")[[1]])
nsweeps <- as.integer(get_arg("nsweeps", 10000))
burnin  <- as.integer(get_arg("burnin", 1000))
thin    <- as.integer(get_arg("thin", 10))
output  <- get_arg("output", "ising_results.csv")

# Make snapshot directory
snapshot_dir <- "grid_snapshot"
if (!dir.exists(snapshot_dir)) dir.create(snapshot_dir)

# Periodic boundary helper
iwrap <- function(i) ((i - 1) %% N) + 1

# Energy function
energy_func <- function(S) {
  E <- 0
  for (i in 1:N) for (j in 1:N) {
    nb_sum <- S[iwrap(i+1), j] + S[iwrap(i-1), j] + S[i, iwrap(j+1)] + S[i, iwrap(j-1)]
    E <- E - S[i,j] * nb_sum / 2
  }
  E
}

results <- list()

# Monte Carlo simulation with progress prints
for (T in temps) {
  cat(sprintf("\n=== Temperature T = %.3f ===\n", T))
  S <- matrix(sample(c(-1,1), N*N, replace=TRUE), nrow=N)
  record_count <- 0
  next_report <- nsweeps / 10

  for (sweep in 1:nsweeps) {
    # Metropolis sweep
    for (i in 1:N) for (j in 1:N) {
      h <- S[iwrap(i+1), j] + S[iwrap(i-1), j] + S[i, iwrap(j+1)] + S[i, iwrap(j-1)]
      dE <- 2 * S[i,j] * h
      if (dE <= 0 || runif(1) < exp(-dE / T)) S[i,j] <- -S[i,j]
    }
    # Progress report every 10%
    if (sweep == floor(next_report)) {
      cat(sprintf("  Completed %d%% (%d/%d sweeps)\n", 
                  round(sweep / nsweeps * 100), sweep, nsweeps))
      next_report <- next_report + nsweeps / 10
    }
    # Recording after burn-in
    if (sweep > burnin && ((sweep - burnin) %% thin == 0)) {
      record_count <- record_count + 1
      M <- mean(S)
      E <- energy_func(S)
      results[[length(results)+1]] <- c(T=T, sweep=sweep, M=M, E=E)

      # write a grid snapshot
      snap_file <- file.path(
        snapshot_dir,
        sprintf("snapshot_%d_T%.3f_sweep%d.csv",
                record_count, T, sweep)
      )
      write.table(
        S,
        file       = snap_file,
        sep        = ",",
        row.names  = FALSE,
        col.names  = FALSE,
        quote      = FALSE
      )
    }
  }
  cat(sprintf("Finished T=%.3f: recorded %d samples.\n", T, record_count))
}

# Combine and write output
df <- do.call(rbind, results)
colnames(df) <- c("T","sweep","M","E")
df <- as.data.frame(df)
write.csv(df, file=output, row.names=FALSE)


cat(sprintf(
  "\nAll simulations complete.\n • Stats in '%s' (%d rows)\n • Snapshots in '%s/' (%d files)\n",
  output, nrow(df),
  snapshot_dir, length(list.files(snapshot_dir))
))
