#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript gwas_gcta_plots.R <gcta_mlma_results.mlma>")
}
infile <- args[1]

cat("Lecture du fichier GWAS :", infile, "\n")
dat <- read.table(infile, header = TRUE, stringsAsFactors = FALSE)

names(dat) <- gsub("\\.", "_", names(dat))
names_low <- tolower(names(dat))

chr_idx <- which(names_low %in% c("chr", "chrom", "chromosome"))
pos_idx <- which(names_low %in% c("bp", "pos", "position", "basepair"))
p_idx   <- which(names_low %in% c("p", "pval", "pvalue", "p_value", "p_wald", "p_lrt"))

if (length(chr_idx) == 0 || length(pos_idx) == 0 || length(p_idx) == 0) {
  stop("Impossible de trouver les colonnes CHR/BP/P. Vérifie le .mlma.")
}

CHR <- as.numeric(dat[[chr_idx[1]]])
BP  <- as.numeric(dat[[pos_idx[1]]])
P   <- as.numeric(dat[[p_idx[1]]])

valid <- !is.na(CHR) & !is.na(BP) & !is.na(P) & P > 0 & P <= 1
CHR <- CHR[valid]
BP  <- BP[valid]
P   <- P[valid]

gwas <- data.frame(CHR = CHR, BP = BP, P = P)
gwas <- gwas[order(gwas$CHR, gwas$BP), ]

# Lambda GC
chisq <- qchisq(1 - gwas$P, df = 1)
lambda_gc <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)

outdir <- "qc_results/gcta"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Lambda GC =", lambda_gc, "\n")
writeLines(sprintf("Lambda_GC\t%0.4f", lambda_gc),
           file.path(outdir, "gcta_lambda_gc.txt"))

# QQ plot
n <- length(gwas$P)
expected <- -log10(ppoints(n))
observed <- -log10(sort(gwas$P))

png(file.path(outdir, "gcta_qqplot.png"),
    width = 960, height = 768)
par(cex.lab = 1.4, cex.main = 1.5)
plot(expected, observed,
     xlab = expression(paste("Attendu -log"[10], "(P)")),
     ylab = expression(paste("Observé -log"[10], "(P)")),
     main = "QQ plot GCTA-MLMA",
     pch = 20)
abline(0, 1, col = "red", lwd = 2)
legend("topleft",
       legend = sprintf("Lambda GC = %.3f", lambda_gc),
       bty = "n")
dev.off()

# Manhattan plot
gwas$CHR <- as.factor(gwas$CHR)
chr_levels <- sort(unique(as.numeric(as.character(gwas$CHR))))
gwas$CHR <- factor(gwas$CHR, levels = chr_levels)

chr_lengths <- tapply(gwas$BP, gwas$CHR, max)
offsets <- c(0, cumsum(as.numeric(chr_lengths))[-length(chr_lengths)])
names(offsets) <- chr_levels

gwas$pos_cum <- gwas$BP + offsets[as.character(gwas$CHR)]
axis_df <- aggregate(pos_cum ~ CHR, data = gwas, FUN = mean)

genomewide <- 5e-8

png(file.path(outdir, "gcta_manhattan.png"),
    width = 1920, height = 1080)
par(mar = c(5, 5, 4, 2))
plot(gwas$pos_cum, -log10(gwas$P),
     pch = 20,
     xlab = "Chromosome",
     ylab = expression(paste("-log"[10], "(P)")),
     main = "Manhattan plot GCTA-MLMA",
     xaxt = "n")
axis(1, at = axis_df$pos_cum, labels = axis_df$CHR)
abline(h = -log10(genomewide), col = "red", lty = 2)
dev.off()

cat("Plots générés :\n")
cat("  -", file.path(outdir, "gcta_manhattan.png"), "\n")
cat("  -", file.path(outdir, "gcta_qqplot.png"), "\n")
cat("  -", file.path(outdir, "gcta_lambda_gc.txt"), "\n")

