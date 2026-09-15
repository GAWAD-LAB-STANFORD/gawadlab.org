# Pack the Sherlock export into the one binary file the browser loads.
# ---------------------------------------------------------------------------
# Inputs (see build.sh for how they are produced):
#   cereb_cells.csv        40,253 cells: t-SNE coordinates, cell type, timepoint, bin
#   cereb_bins.csv         996 grid squares over the embedding
#   cereb_expr.bin.gz      dense uint8, bins x genes, column-major
#   cereb_expr_genes.txt   column names for that matrix
#   cereb_expr_meta.csv    n_bins, n_genes, and the value the uint8 range maps to
#   cereb_groups.csv.gz    per-gene mean and detection rate per cell type / timepoint
#   cereb_group_sizes.csv  how many cells each of those groups holds
#
# Text parsing is what makes a WebAssembly app sit on a spinner, so everything
# lands in one .rds that R deserialises. Group statistics are stored as scaled
# integers rather than doubles - same precision as the source CSV, far smaller.
#
#   Rscript prep_bundle.R DATA_DIR OUT.rds
# ---------------------------------------------------------------------------

build_bundle <- function(D) {
  rp <- function(f) file.path(D, f)

  meta  <- read.csv(rp("cereb_expr_meta.csv"))
  genes <- readLines(rp("cereb_expr_genes.txt"))
  stopifnot(length(genes) == meta$n_genes)

  con  <- gzfile(rp("cereb_expr.bin.gz"), "rb")
  bits <- readBin(con, "raw", n = meta$n_bins * meta$n_genes); close(con)
  stopifnot(length(bits) == meta$n_bins * meta$n_genes)
  gmat <- matrix(bits, nrow = meta$n_bins, dimnames = list(NULL, genes))

  cells <- read.csv(rp("cereb_cells.csv"), stringsAsFactors = FALSE)
  bins  <- read.csv(rp("cereb_bins.csv"),  stringsAsFactors = FALSE)
  gsz   <- read.csv(rp("cereb_group_sizes.csv"), stringsAsFactors = FALSE)
  grp   <- read.csv(gzfile(rp("cereb_groups.csv.gz")), check.names = FALSE,
                    stringsAsFactors = FALSE)

  CT <- sub("^ct_mean\\|", "", grep("^ct_mean\\|", names(grp), value = TRUE))
  TP <- sub("^tp_mean\\|", "", grep("^tp_mean\\|", names(grp), value = TRUE))
  # E10..E18 then P0..P10 - sorted as text, P10 would land between P0 and P4
  TP <- TP[order(match(TP, c(paste0("E", 8:22), paste0("P", 0:30))), TP)]

  as_int <- function(pre, keys, scale) {
    m <- as.matrix(grp[, paste0(pre, keys), drop = FALSE])
    storage.mode(m) <- "double"
    m <- matrix(as.integer(round(m * scale)), nrow = nrow(m),
                dimnames = list(grp$gene, keys))
    m
  }
  CT_MEAN <- as_int("ct_mean|", CT, 1e4); CT_PCT <- as_int("ct_pct|", CT, 1e2)
  TP_MEAN <- as_int("tp_mean|", TP, 1e4); TP_PCT <- as_int("tp_pct|", TP, 1e2)

  # per-gene summary, so the browser never recomputes it
  ord  <- t(apply(CT_MEAN, 1, sort, decreasing = TRUE))
  summ <- data.frame(
    gene        = grp$gene,
    peak_type   = CT[max.col(CT_MEAN, "first")],
    specificity = round(ord[, 1] / pmax(rowSums(CT_MEAN), 1), 3),
    fold_2nd    = round(ord[, 1] / pmax(ord[, 2], 1), 2),
    peak_time   = TP[max.col(TP_MEAN, "first")],
    max_pct     = round(apply(CT_PCT, 1, max) / 1e2, 1),
    mapped      = grp$gene %in% genes,
    stringsAsFactors = FALSE)

  # the cell type that dominates each grid square, for the comparison scatter
  tb <- table(cells$bin, cells$cell_type)
  bins$cell_type <- colnames(tb)[max.col(tb, "first")][
    match(bins$bin, as.integer(rownames(tb)))]

  list(cells = data.frame(D1 = cells$D1, D2 = cells$D2,
                          cell_type = factor(cells$cell_type, levels = CT),
                          timepoint = factor(cells$timepoint, levels = TP),
                          bin = as.integer(cells$bin)),
       bins = bins, gmat = gmat, genes = genes, vmax = meta$vmax,
       CT = CT, TP = TP, CT_MEAN = CT_MEAN, CT_PCT = CT_PCT,
       TP_MEAN = TP_MEAN, TP_PCT = TP_PCT, MSCALE = 1e4, PSCALE = 1e2,
       summ = summ, gsz = gsz)
}

if (!interactive() && sys.nframe() == 0L) {
  a <- commandArgs(TRUE)
  B <- build_bundle(a[1])
  saveRDS(B, a[2], compress = "gzip")
  cat(sprintf("bundle -> %s (%.1f MB) | %s cells | %s bins | %s genes mapped, %s with statistics\n",
              a[2], file.size(a[2]) / 1e6, format(nrow(B$cells), big.mark = ","),
              format(nrow(B$bins), big.mark = ","), format(length(B$genes), big.mark = ","),
              format(nrow(B$summ), big.mark = ",")))
}
