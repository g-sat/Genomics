# Partitioning versus cost per sample.
# Reads the table written by experiments/cost_per_sample.py and redraws the
# paper figure. The CSV is the simulator comparison of the nuclear-family
# oracle against the all-samples oracle. HDFS bytes per sample are the
# one-host select_05 measurement.
#
#   Rscript analysis/partition_cost.R

args_root <- if (file.exists("results/cost_per_sample/ratios.csv")) {
  "."
} else if (file.exists("../results/cost_per_sample/ratios.csv")) {
  ".."
} else {
  stop("run python experiments/cost_per_sample.py first")
}

ratios <- read.csv(file.path(args_root, "results/cost_per_sample/ratios.csv"))
hdfs <- read.csv(file.path(args_root, "results/cost_per_sample/hdfs_per_sample.csv"))

cohorts <- unique(ratios$cohort)
storage <- sapply(cohorts, function(name) median(ratios$storage_ratio[ratios$cohort == name]))
query <- sapply(cohorts, function(name) median(ratios$query_ratio[ratios$cohort == name]))
shuffle_kb <- sapply(cohorts, function(name) {
  median(ratios$shuffle_bytes_per_sample[ratios$cohort == name]) / 1000
})

out <- file.path(args_root, "results/cost_per_sample/cost_per_sample_r.png")
png(out, width = 1000, height = 430)
par(mfrow = c(1, 2), mar = c(6, 5, 3, 1))
bar_x <- barplot(
  rbind(storage, query),
  beside = TRUE,
  names.arg = cohorts,
  col = c("#4c78a8", "#f58518"),
  ylim = c(0, 1.7),
  ylab = "Nuclear family / all-samples",
  main = "Cost per sample of each partition",
  las = 2
)
abline(h = 1, lwd = 1)
legend("topright", c("Storage seconds / sample", "Fit-query seconds / sample"),
       fill = c("#4c78a8", "#f58518"), bty = "n", cex = 0.8)
barplot(
  shuffle_kb,
  names.arg = cohorts,
  col = "#4c78a8",
  ylab = "Shuffle on the fit mix (KB / sample)",
  main = "Hereditary check that crosses the partition",
  las = 2
)
dev.off()

message("wrote ", out)
message("one-host HDFS, held-out hereditary-check bytes per sample:")
print(hdfs[c("seed", "partition", "check_bytes_per_sample", "primary_per_sample")])
