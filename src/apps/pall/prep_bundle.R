# Pack the pALL manuscript's figure-level source data into one binary bundle.
# ---------------------------------------------------------------------------
# Reads 08_Source_Data from the Nature Genetics submission package and writes the
# single .rds the browser loads. Everything here is plot-level data that already
# ships with the paper; no raw sequence, no patient-identifiable fields.
#
#   Rscript prep_bundle.R SOURCE_DATA_DIR out.rds
# ---------------------------------------------------------------------------

build_bundle <- function(D) {
  rp  <- function(...) file.path(D, ...)
  rc  <- function(...) read.csv(rp(...), stringsAsFactors = FALSE, check.names = FALSE,
                                fileEncoding = "UTF-8-BOM")
  rt  <- function(...) read.delim(rp(...), stringsAsFactors = FALSE, check.names = FALSE)

  # --- 1. bulk vs single-cell mutation burden (Fig 3D) --------------------
  burden <- rc("Figure_3", "Somatic_DP5_cells.csv")
  names(burden) <- make.names(names(burden), unique = TRUE)
  burden <- burden[, c("Sample", "Color", "Total.Somatic", "Unique", "Shared",
                       "Som..Mb", "Sensitivity", "Cor.Som.Mb", "Cor..Total")]
  names(burden) <- c("assay", "patient", "total_somatic", "unique", "shared",
                     "som_per_mb", "sensitivity", "corrected_som_per_mb", "corrected_total")
  burden$patient <- sub("^ID", "", burden$patient)
  burden <- burden[burden$assay %in% c("Bulk", "Single"), ]

  # --- 2. where pALL sits among paediatric tumours (Fig 3E reference) -----
  pan <- rc("Figure_3", "SNV_per_Mb.csv")
  names(pan) <- make.names(names(pan), unique = TRUE)
  pan <- data.frame(sample = pan$Sample, cancer_type = pan$Cancer.Type,
                    snvs = pan$SNVs, snv_per_mb = pan$SNVs_Per_Mb,
                    total_per_mb = pan$Total_Mutations_Per_Mb, stringsAsFactors = FALSE)
  pan <- pan[is.finite(pan$snv_per_mb), ]

  # --- 3. activating RAS mutations found by error-corrected sequencing ----
  ras <- rt("Figure_1", "Ras_AF2.tsv")
  ras$AF <- as.numeric(ras$AF)
  ras <- ras[is.finite(ras$AF), ]

  # --- 4. ex vivo drug response, five patients (Fig 6C) -------------------
  #     rows are "<patient> | <patient>_<condition>_<replicate>", columns are
  #     "<GENE> p.<change>", values are mutant allele frequency in percent
  dm <- rt("Figure_6", "Tamara_Fig6B", "2026-05-28_ALL_Vast_fixed.clustered_matrix.tsv")
  lab <- sub("^.*\\|\\s*", "", dm$sample)
  drug <- data.frame(
    patient   = sub("\\s*\\|.*$", "", dm$sample),
    condition = sub("^[0-9]+_", "", sub("_[0-9]+$", "", lab)),
    replicate = ifelse(grepl("_[0-9]+$", lab), sub("^.*_", "", lab), ""),
    stringsAsFactors = FALSE)
  dmat <- as.matrix(dm[, setdiff(names(dm), "sample"), drop = FALSE])
  storage.mode(dmat) <- "double"
  keep <- colSums(dmat > 0, na.rm = TRUE) > 0          # drop all-zero mutations
  dmat <- dmat[, keep, drop = FALSE]

  # --- 5. ex vivo drug response, SJETV077 across nine conditions (Fig 6B) -
  sj <- rt("Figure_6", "all_mutation_mean.tsv")
  names(sj)[names(sj) == "x"] <- "af"
  sj$Treatment <- sub("^annotated\\.LTIVC-SJETV077-TB-05-2296-", "", sj$Treatment)
  sj$gene <- sub(":.*$", "", sj$Mutation)
  sj$change <- ifelse(grepl("p\\.", sj$Mutation), sub("^.*p\\.", "p.", sj$Mutation), "")

  # --- 6. single cells before and after treatment (Fig 7A) ----------------
  cells <- rc("Figure_7", "CNV_CSF_Plot.csv")
  clone <- rt("Figure_7", "Tamara_4295_source_data", "cnv", "4272-4295-CNV_Cluser_ID.tsv")
  clone$sample <- sub("[.].*$", "", clone$GeneID)
  clone$cell   <- sub("^[^.]*[.]", "", clone$GeneID)
  cells$clone  <- clone$Cluster[match(paste(cells$TimeLine, cells$Cell_ID),
                                      paste(clone$sample, clone$cell))]
  # 4272 holds no cells from clusters 6 or 7, which is what makes those emergent
  em <- setdiff(unique(clone$Cluster[clone$sample == "4295"]),
                unique(clone$Cluster[clone$sample == "4272"]))
  cells$timepoint <- factor(ifelse(cells$TimeLine == "4272", "Pretreatment", "Post-treatment"),
                            levels = c("Pretreatment", "Post-treatment"))
  geno <- rc("Figure_7", "_B.csv")
  gmat <- as.matrix(geno[, -1, drop = FALSE]); storage.mode(gmat) <- "integer"
  rownames(gmat) <- geno$SampleName

  # --- 7. gene-level recurrence and predicted pathogenicity (Fig 7B, 7C) --
  rec <- rc("Figure_7", "panels_B_C", "panel_2axis.csv")
  am  <- rc("Figure_7", "panels_B_C", "observed_variants_am2.csv")
  am$am_pathogenicity <- suppressWarnings(as.numeric(am$am_pathogenicity))
  am <- am[is.finite(am$am_pathogenicity), ]
  # the published panels drop three genes inadvertently included in the control set
  am  <- am[!(am$gene %in% c("JAK2", "FLT3", "CDKN2A")), ]
  rec <- rec[!(rec$gene %in% c("JAK2", "FLT3", "CDKN2A")), ]

  # --- 9. allele frequency before and after treatment (Fig 7 paired samples) -
  #     Each variant appears twice in the ConDoR read-count files, once suffixed
  #     with its gene and once with _NA. The _NA column is an empty placeholder -
  #     zero alt reads and no total depth anywhere - so only the gene-named
  #     columns carry data. Two multi-allelic sites have no gene column at all.
  ar <- rc("Figure_7", "ConDoR_alt_readscount.csv")
  tr <- rc("Figure_7", "ConDoR_total_readscount.csv")
  stopifnot(identical(names(ar), names(tr)), identical(ar[[1]], tr[[1]]))
  vc <- names(ar)[-1]; vc <- vc[sub("^.*_", "", vc) != "NA"]
  A <- as.matrix(ar[, vc, drop = FALSE]); TT <- as.matrix(tr[, vc, drop = FALSE])
  storage.mode(A) <- "double"; storage.mode(TT) <- "double"
  stopifnot(all(colSums(TT, na.rm = TRUE) > 0))
  tp <- sub("-.*", "", ar[[1]])
  pseudobulk <- function(sel) {
    al <- colSums(A[sel, , drop = FALSE], na.rm = TRUE)
    to <- colSums(TT[sel, , drop = FALSE], na.rm = TRUE)
    ifelse(to > 0, al / to, NA_real_)
  }
  ncov <- function(sel) colSums(!is.na(TT[sel, , drop = FALSE]) & TT[sel, , drop = FALSE] > 0)
  i1 <- tp == "4272"; i2 <- tp == "4295"
  af <- data.frame(
    variant    = vc,
    gene       = sub("^.*_", "", vc),
    locus      = sub("_[^_]*$", "", vc),
    pre        = pseudobulk(i1),  post       = pseudobulk(i2),
    cells_pre  = ncov(i1),        cells_post = ncov(i2),
    depth_pre  = colSums(TT[i1, , drop = FALSE], na.rm = TRUE),
    depth_post = colSums(TT[i2, , drop = FALSE], na.rm = TRUE),
    stringsAsFactors = FALSE)
  af$delta <- af$post - af$pre
  af <- af[is.finite(af$pre) & is.finite(af$post), ]
  rownames(af) <- NULL

  # --- 8. phylogenies (Fig 5E/5G maximum likelihood, Fig 7A clone tree) ---
  #     Newick is kept as text and parsed in the browser, so the bundle carries
  #     kilobytes rather than a serialised tree object.
  PH <- c("4295", "445", "417", "4084", "Invitro")
  trees <- list()
  for (q in PH) {
    f <- rp("shared_phycall_Figures_4_5_S6_S7",
            sprintf("CellPhy.%s.GT10+FO+E.nobulks.noCNVs.raxml.support", q))
    if (file.exists(f)) trees[[q]] <- paste(readLines(f, warn = FALSE), collapse = "")
  }
  trees[["clone"]] <- paste(readLines(rp("Figure_7", "_tree.newick"), warn = FALSE), collapse = "")

  # one annotation row per cell, keyed on a normalised id, because the trees
  # write 4295_A10 / 4295.F1 while the tables write 4295-A10
  norm <- function(x) gsub("[._]", "-", x)
  tips <- data.frame(
    cell      = norm(cells$Index),
    sample    = as.character(cells$TimeLine),
    timepoint = as.character(cells$timepoint),
    clone     = ifelse(is.na(cells$clone), "unassigned", cells$clone),
    chr4      = cells$Chr4_Deletion,
    chr6      = cells$Chr6_1_Deletion,
    stringsAsFactors = FALSE)

  list(trees = trees, tips = tips, af = af,
       burden = burden, pan = pan, ras = ras,
       drug = drug, dmat = dmat, sj = sj,
       cells = cells, gmat = gmat, emergent = em,
       rec = rec, am = am)
}

if (!interactive() && sys.nframe() == 0L) {
  a <- commandArgs(TRUE)
  B <- build_bundle(a[1])
  saveRDS(B, a[2], compress = "gzip")
  cat(sprintf("bundle -> %s (%.2f MB)\n", a[2], file.size(a[2]) / 1e6))
  cat(sprintf("  burden      %d rows (%d patients, bulk vs single)\n", nrow(B$burden),
              length(unique(B$burden$patient))))
  cat(sprintf("  paediatric  %d tumours across %d cancer types\n", nrow(B$pan),
              length(unique(B$pan$cancer_type))))
  cat(sprintf("  RAS         %d mutations in %d patients\n", nrow(B$ras),
              length(unique(B$ras$Patient))))
  cat(sprintf("  drug        %d samples x %d mutations, %d patients, conditions: %s\n",
              nrow(B$drug), ncol(B$dmat), length(unique(B$drug$patient)),
              paste(sort(unique(B$drug$condition)), collapse = ", ")))
  cat(sprintf("  SJETV077    %d mutations x %d conditions\n",
              length(unique(B$sj$Mutation)), length(unique(B$sj$Treatment))))
  cat(sprintf("  cells       %d (%s), %d clones, emergent: %s\n", nrow(B$cells),
              paste(names(table(B$cells$timepoint)), table(B$cells$timepoint),
                    sep = "=", collapse = ", "),
              length(unique(na.omit(B$cells$clone))), paste(B$emergent, collapse = ", ")))
  cat(sprintf("  genotypes   %d cells x %d variants\n", nrow(B$gmat), ncol(B$gmat)))
  cat(sprintf("  genes       %d with recurrence, %d scored missense variants\n",
              nrow(B$rec), nrow(B$am)))
  cat("  trees       ")
  for (n in names(B$trees)) cat(sprintf("%s(%d chars) ", n, nchar(B$trees[[n]])))
  cat(sprintf("\n  tip labels  %d annotated cells\n", nrow(B$tips)))
  cat(sprintf("  before/after %d variants with pseudobulk VAF at both timepoints; largest rise %s %+.3f, largest fall %s %+.3f\n",
              nrow(B$af), B$af$gene[which.max(B$af$delta)], max(B$af$delta),
              B$af$gene[which.min(B$af$delta)], min(B$af$delta)))
}
