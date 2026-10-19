# Exports qqman's example data and renders reference plots with real R/qqman,
# for side-by-side comparison with qqman-rs output.
#
#   Rscript validation/reference.R     # -> validation/data/*, validation/reference/*.{pdf,png}
#
# Paths are relative to this script's folder, so it can be run from anywhere.
suppressPackageStartupMessages(library(qqman))
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
dir.create(file.path(here, "data"), showWarnings = FALSE)
dir.create(file.path(here, "reference"), showWarnings = FALSE)
write.table(gwasResults, file.path(here, "data", "gwasResults.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

snpsOfInterest <- snpsOfInterest
writeLines(snpsOfInterest, file.path(here, "data", "snpsOfInterest.txt"))

plots <- list(
  manhattan = function() manhattan(gwasResults),
  manhattan_main = function() manhattan(gwasResults, main = "Manhattan Plot", ylim = c(0, 10),
                                        cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"),
                                        suggestiveline = FALSE, genomewideline = FALSE),
  manhattan_chr3 = function() manhattan(subset(gwasResults, CHR == 3)),
  manhattan_highlight = function() manhattan(gwasResults, highlight = snpsOfInterest),
  manhattan_annotate = function() manhattan(gwasResults, annotatePval = 0.01),
  qq = function() qq(gwasResults$P)
)
for (name in names(plots)) {
  # qqman-rs's default sizes: manhattan 14 x 7 in, qq 7 x 7 in; PNGs at 150 dpi.
  width <- if (startsWith(name, "manhattan")) 14 else 7
  pdf(file.path(here, "reference", paste0(name, ".pdf")), width = width, height = 7)
  plots[[name]]()
  invisible(dev.off())
  png(file.path(here, "reference", paste0(name, ".png")), width = width, height = 7, units = "in", res = 150)
  plots[[name]]()
  invisible(dev.off())
}
# R's PNG output depends on the device (quartz on macOS, cairo elsewhere), so
# record which one made these references.
writeLines(c(R.version.string, paste("png() type:", getOption("bitmapType")),
             paste("qqman", packageVersion("qqman"))),
           file.path(here, "reference", "png-device.txt"))
