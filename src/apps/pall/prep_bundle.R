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

  # --- 10. per-cell metadata for every patient -----------------------------
  #     The Figure 7 tables cover one patient's paired samples only. These two
  #     sources cover all four patients plus the in vitro benchmark, so the
  #     phylogenies of 417, 445 and 4084 can be coloured by something real
  #     instead of rendering uniformly grey.
  PP <- file.path(D, "shared_phycall_Figures_4_5_S6_S7")
  # stored as a data.frame; as.vector() on one of those returns columns rather
  # than values, so coerce to a numeric matrix once, here
  sig <- as.matrix(readRDS(file.path(PP, "Signatures.rds")))  # 191 x 30 exposures
  storage.mode(sig) <- "double"
  nsnv <- as.data.frame(readRDS(file.path(PP, "Signatures.suppl.info.rds")))
  sig <- sig[!grepl("-(germline|internal)$", rownames(sig)), , drop = FALSE]
  top_sig <- colnames(sig)[max.col(sig, "first")]
  top_sig[rowSums(sig) == 0] <- NA
  sigdf <- data.frame(cell = gsub("[._]", "-", rownames(sig)),
                      patient = sub("[-_].*", "", rownames(sig)),
                      top_signature = sub("Signature\\.", "SBS", top_sig),
                      top_share = round(apply(sig, 1, max), 3),
                      stringsAsFactors = FALSE)
  sigdf$n_snv <- nsnv$NumSNVs[match(rownames(sig), nsnv$Sample)]

  ado <- do.call(rbind, lapply(c("4295", "4084", "417", "445"), function(q) {
    f <- file.path(PP, sprintf("ADO.%s.rds", q))
    if (!file.exists(f)) return(NULL)
    a <- readRDS(f)
    data.frame(cell = gsub("[._]", "-", paste0(q, "-", trimws(a$Sample))),
               patient = q, ado = suppressWarnings(as.numeric(a$ado)),
               depth = suppressWarnings(as.numeric(a$seq)),
               method = a$method, stringsAsFactors = FALSE)
  }))
  ado <- ado[!grepl("bulk", ado$cell, ignore.case = TRUE), ]

  cellmeta <- merge(sigdf, ado[, c("cell", "ado", "depth", "method")], by = "cell", all = TRUE)
  cellmeta$patient[is.na(cellmeta$patient)] <- sub("[-_].*", "", cellmeta$cell[is.na(cellmeta$patient)])
  # the exposure matrix itself, for the per-patient signature panel
  rownames(sig) <- gsub("[._]", "-", rownames(sig))
  colnames(sig) <- sub("Signature\\.", "SBS", colnames(sig))
  sig <- sig[, colSums(sig) > 0, drop = FALSE]
  stopifnot(is.matrix(sig), is.numeric(sig))

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

  # Per-cell mutation composition for the exome tree: each cell's alternate reads
  # split across the panel's genes. Signature fits exist only for the
  # post-treatment cells, so they cannot colour a before-and-after tree; these
  # read counts cover 110 of its 113 tips including 29 of 30 pretreatment cells.
  gene_of <- sub("^.*_", "", vc)
  mut <- t(rowsum(t(A), gene_of, na.rm = TRUE))          # cells x gene, alt reads
  rownames(mut) <- gsub("[._]", "-", ar[[1]])
  mut <- mut[, colSums(mut, na.rm = TRUE) > 0, drop = FALSE]
  storage.mode(mut) <- "double"

  # --- 11. induction therapy across all four patients ----------------------
  #     Two different strategies built the two groups of phylogenies, and they
  #     carry "before and after" in two different places.
  #
  #     (a) Figure 5: PTA single-cell WHOLE GENOMES, trees built by CellPhy.
  #         Four patients, 165 cells. Every one of those cells was drawn AFTER
  #         four weeks of induction therapy - the manuscript's own words are
  #         "four high-risk pALL patients who had undergone four weeks of
  #         induction therapy ... but had detectable MRD cells". So the trees
  #         themselves have no before/after axis.
  #     (b) The before/after for those four patients is BULK, in
  #         Mutect<patient>.filtered.rds: every exome-mapped tree mutation
  #         re-measured in the diagnostic bulk and again in the remission bulk.
  #     (c) Figure 7: single-cell EXOME, tree built by ConDoR. One patient only
  #         (4295), and there the before/after is per cell. That is the tree the
  #         rest of this tab draws.
  #
  #     Sample naming is not documented in a table, so it is derived: strip the
  #     repeated "<patient>-" prefix, and what remains is "<accession><suffix>".
  #     A cell leaves no accession behind (417-A10 -> "A10"), a bulk does
  #     (417-368B -> accession 368, suffix B). The notebook's own rule at
  #     Figures.BALL.PTA.Rmd:1508 is timing = if the accession equals the patient
  #     id then "remission" else "diagnosis"; the patient is named for its
  #     remission accession. The VAFs below confirm it independently in all four
  #     patients: the other accession carries these mutations at ~0.3-0.5 VAF
  #     (clonal disease) and the patient-named one at ~0.01-0.03 (residual).
  #     Suffix NonB is the sorted non-blast fraction, i.e. a normal control, so
  #     it is kept separately and never averaged into the tumour columns.
  wgs <- list(); bulk <- list()
  for (q in c("417", "445", "4084", "4295")) {
    td <- readRDS(rp("shared_phycall_Figures_4_5_S6_S7",
                     sprintf("TreeMutWithZeros.%s.genome.exomemapped.rds", q)))
    ph <- td@phylo
    dd <- as.data.frame(td@data)
    nt <- ape::Ntip(ph)

    # Identify each annotated branch by the SET OF TIPS beneath it rather than by
    # its node index, so the annotation survives the newick round trip into webR
    # (ape renumbers nodes from the edge matrix, not from the file).
    # descendant tips, by walking the edge matrix - avoids a phangorn dependency
    kids <- split(ph$edge[, 2], ph$edge[, 1])
    tips_under <- function(n) {
      if (n <= nt) return(n)
      out <- integer(0); stack <- as.integer(kids[[as.character(n)]])
      while (length(stack)) {
        v <- stack[1]; stack <- stack[-1]
        if (v <= nt) out <- c(out, v) else stack <- c(stack, as.integer(kids[[as.character(v)]]))
      }
      sort(out)
    }
    desc <- function(n) paste(sort(ph$tip.label[tips_under(n)]), collapse = "|")
    ann <- dd[!is.na(dd$mutList) & nzchar(dd$mutList), c("node", "mutList")]
    ann <- ann[order(ann$node), ]
    wgs[[q]] <- list(
      nwk = ape::write.tree(ph),
      mut = if (nrow(ann)) data.frame(
        key   = vapply(ann$node, desc, character(1)),
        n     = vapply(ann$node, function(n) length(tips_under(n)), integer(1)),
        label = gsub("\n", "; ", ann$mutList),
        stringsAsFactors = FALSE) else
        data.frame(key = character(0), n = integer(0), label = character(0)))

    m <- readRDS(rp("shared_phycall_Figures_4_5_S6_S7",
                    sprintf("Mutect%s.filtered.rds", q)))
    rest <- as.character(m$sample)
    repeat { r2 <- sub(paste0("^", q, "-"), "", rest); if (identical(r2, rest)) break; rest <- r2 }
    acc <- sub("[^0-9].*$", "", rest)          # empty for a single cell
    suf <- sub("^[0-9]+", "", rest)
    isbulk <- nzchar(acc)
    stopifnot(any(isbulk), any(!isbulk))
    timing <- ifelse(acc == q, "after", "before")
    control <- grepl("NonB", suf)
    agg <- function(keep) {
      z <- m[isbulk & keep, ]
      if (!nrow(z)) return(NULL)
      # one patient has two diagnostic aliquots (T1/T2); pool their reads rather
      # than averaging two VAFs computed at different depths
      data.frame(id = z$id, alt = z$VAF * z$DP, dp = z$DP, stringsAsFactors = FALSE)
    }
    pool <- function(z) if (is.null(z)) NULL else {
      a <- rowsum(z[, c("alt", "dp")], z$id)
      data.frame(id = rownames(a), vaf = ifelse(a$dp > 0, a$alt / a$dp, NA_real_),
                 dp = a$dp, stringsAsFactors = FALSE)
    }
    bf <- pool(agg(timing == "before" & !control))
    af2 <- pool(agg(timing == "after"  & !control))
    nb <- pool(agg(control))
    stopifnot(!is.null(bf), !is.null(af2))
    ids <- sort(unique(c(bf$id, af2$id)))
    ann1 <- m[match(ids, m$id), c("id", "Gene.refGene", "ExonicFunc.refGene", "CHROM:POS")]
    bulk[[q]] <- data.frame(
      patient  = q,
      id       = ids,
      gene     = as.character(ann1$Gene.refGene),
      change   = sub("^.*:p\\.", "", ids),
      effect   = sub("_SNV$", "", as.character(ann1$ExonicFunc.refGene)),
      locus    = as.character(ann1[["CHROM:POS"]]),
      before   = bf$vaf[match(ids, bf$id)],
      after    = af2$vaf[match(ids, af2$id)],
      dp_before= bf$dp[match(ids, bf$id)],
      dp_after = af2$dp[match(ids, af2$id)],
      normal   = if (is.null(nb)) NA_real_ else nb$vaf[match(ids, nb$id)],
      stringsAsFactors = FALSE)
  }
  bulk <- do.call(rbind, bulk); rownames(bulk) <- NULL
  bulk$change[!grepl(":p\\.", bulk$id)] <- ""
  stopifnot(nrow(bulk) > 0, all(is.finite(bulk$before) | is.finite(bulk$after)))

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

  list(trees = trees, tips = tips, af = af, cellmeta = cellmeta, sig = sig, mut = mut,
       wgs = wgs, bulk = bulk,
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
  cat(sprintf("  cell meta   %d cells across %s\n", nrow(B$cellmeta),
              paste(sort(unique(B$cellmeta$patient)), collapse = ", ")))
  cat(sprintf("  signatures  %d cells x %d active COSMIC signatures; most common top signature %s\n",
              nrow(B$sig), ncol(B$sig),
              names(sort(table(B$cellmeta$top_signature), decreasing = TRUE))[1]))
  cat(sprintf("  mutation pie %d cells x %d genes from alternate reads; %d cells carry any\n",
              nrow(B$mut), ncol(B$mut), sum(rowSums(B$mut, na.rm = TRUE) > 0)))
  cat(sprintf("  before/after %d variants with pseudobulk VAF at both timepoints; largest rise %s %+.3f, largest fall %s %+.3f\n",
              nrow(B$af), B$af$gene[which.max(B$af$delta)], max(B$af$delta),
              B$af$gene[which.min(B$af$delta)], min(B$af$delta)))
}
