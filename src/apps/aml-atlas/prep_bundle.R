# Build the AML atlas app's data bundle.
# ---------------------------------------------------------------------------
# The app reads one binary bundle. This script builds it from the manuscript
# Source Data, and writes the per-gene expression chunks the app fetches at
# runtime. It replaces the short prep_bundle.R in the submission package, which
# built 19 of the bundle's 36 objects; the other 17 had no build script on disk
# and were produced by scripts that are gone.
#
#   Rscript src/apps/aml-atlas/prep_bundle.R SOURCE_DATA [options]
#
#     --cc DIR        CellChat results: g00.rds .. g74.rds, one per group
#     --expr PATH     all-gene per-bin expression (see "Inputs not in the
#                     package" below). .csv/.csv.gz with bin,gene,mean, or
#                     .rds holding a numeric matrix [bins x genes]
#     --genes-out DIR write the runtime expression chunks here
#     --carry PATH    take all_genes/all_gmax/chunk_size from an existing
#                     bundle instead of rebuilding them from --expr
#     --out PATH      bundle destination (default app/bundle.rds)
#     --check PATH    after building, diff every object against this bundle
#
# WHERE EACH OBJECT COMES FROM
# ---------------------------------------------------------------------------
# All of the following were verified value-by-value against the bundle that
# shipped on 2026-10-21. Rounding is stated where the bundle rounds.
#
# From SOURCE_DATA, unchanged:
#   screen   Figure_5/fig5A_top15.csv            funnel   Figure_5/fig5A_funnel.csv
#   h2h      Figure_5/fig5C_headtohead.csv       ls_expr  Figure_5/fig5D_ls.csv
#   cyto     Figure_5/fig5D_cyto.csv             meta_ped Figure_5/fig5F_meta_all.csv
#   meta_ad  Figure_5/fig5F_adult_beataml_meta.csv
#   sc_ad    Figure_S23/adult_scrna_overall.csv  radar    Figure_6/fig6D_radar.csv
#   tox_hema Figure_6/fig6tox_hema.csv           tox_org  Figure_6/fig6tox_organ.csv
#   tox_ct   Figure_6/fig6tox_celltypes.csv      combos   Figure_7/fig7_allpairs.csv
#   lsspec   Figure_S16/funnel_ls_specific.csv
#   lscomp   Figure_S16/funnel_ls_specific_composition.csv
#   ebins    Figure_2/fig2A_bins.csv
#   tcr      Figure_4/v2/tcr_percell_exh.csv
#   tnk_de   Figure_4/v2/tnk_aml_vs_hbm_de_bySubset.csv
#
# From SOURCE_DATA, derived:
#   emb      Figure_2/fig2A_umap_cells.csv. The file carries no sample column;
#            the barcode suffix after the last "-" is a 0-based sample index,
#            and Figure_4/v2/tnk_subset_percell.csv maps that index 1:1 to a
#            SampleID. marrow follows from the label: 0_Healthy* is a donor.
#   tnk      Figure_S12/figS11_tnk_surface_percell.csv, which carries the
#            embedding, the subset, the marrow and all eight CITE-seq gates.
#            bin is the whole-atlas grid bin, joined from emb on barcode.
#   tnk_markers  the *_pos columns of tnk, in file order.
#   de_all   Figure_2/v2_recluster/group_DE_wholecohort.csv, restricted from
#            27,346 rows to the atlas's own genes; log2FC rounded to 4 dp.
#   ls_up    Figure_3/v3/LS_upgenes_wilcoxon.csv; log2FC rounded to 3 dp.
#            NOTE: the package also ships Figure_2/v2_recluster/
#            LS48_upgenes_wilcoxon.csv, a DIFFERENT version (22,917 rows /
#            8,467 genes against 22,220 / 8,289). The app uses Figure_3/v3.
#   cyto_bulk  Figure_S20/figS_cyto_pediatric_violin.csv and
#            figS_cyto_adult_violin.csv; log2cpm rounded to 4 dp. Four
#            identical copies of the paediatric file ship under Figure_S17,
#            S20, S21 and S22.
#   gmat/gmax/genes_mapped  Figure_2/fig2A_gene_bin_expr.csv, quantised to
#            uint8 against its own maximum. The CSV holds 1,311 genes; the
#            shipped bundle inlined 1,317, the extra six being the clinical
#            comparators ABCA7, CD33, CD38, CLEC12A, FLT3 and IL3RA. Those six
#            are added here only when --expr is given. Leaving them out costs
#            one runtime fetch each, because gene_bins() falls back to the
#            chunks for any gene absent from genes_mapped.
#
# From --cc, the CellChat run output:
#   cc_groups/cc_leuk/cc_normal  the 49 states, then the unassigned leukemic
#            cells, then the non-leukemic cell types in sorted order. This
#            order fixes the gNN.rds filenames the app fetches, so it is
#            constructed rather than read.
#   cc_pairs  pooled from the per-group files and DEDUPLICATED: each file holds
#            every link where its group is source OR target, so pooling all 75
#            double-counts. Dedup on source+target+interaction gives 933,682
#            significant links, the run total.
#   cc_gene_map/cc_db_genes  CellChatDB.human. Not derivable from the link
#            files: the map covers all 3,233 DB interactions where only 2,310
#            occur in the data, and it expands complexes (INHBABB becomes
#            INHBA and INHBB), which splitting the interaction name does not.
#            Written out as cc_gene_map.csv beside the bundle so a later run
#            needs no CellChat install.
#
# INPUTS NOT IN THE PACKAGE
# ---------------------------------------------------------------------------
# Two things cannot be rebuilt from SOURCE_DATA, and this script will not
# invent them. It stops unless you pass them or explicitly opt out.
#
#   1. The all-gene per-bin expression behind --expr. The package deposits
#      per-bin values for the 1,311 screen genes only (fig2A_gene_bin_expr.csv);
#      the atlas quantifies 21,895. The meancpm.npy files are cluster-level
#      (64 x 27,346), not per-bin, so they are not this input. Regenerate it
#      from the count matrix at the 90x90 grid the embedding already defines,
#      or pass --carry to reuse the values from an existing bundle.
#   2. The CellChat run behind --cc. 11h38m on 480 GB over 75 groups; it is a
#      deliverable, not something this script recomputes.
#
# Without --expr or --carry the "Any gene in the atlas" scope is unavailable.
# Without --cc the Cell talk tab is unavailable. Either way the script says so
# and refuses unless you pass --allow-partial.
# ---------------------------------------------------------------------------

say <- function(...) cat(sprintf(...), "\n", sep = "")

# The four leukemic annotations. Everything else in the atlas is a normal
# marrow population, which is what makes cc_normal 25 of the 29 cell types.
LEUK_TYPES <- c("AML", "AML-PCNA", "AML-MKI67", "AML-CD1C")

# uint8 quantisation, the form the browser reads. Returned with the scale so a
# caller cannot lose it: value = as.integer(raw) * gmax / 255.
quantise <- function(M, gmax = max(M)) {
  stopifnot(is.matrix(M), !anyNA(M), gmax > 0)
  Q <- matrix(as.integer(round(M / gmax * 255)), nrow = nrow(M), dimnames = dimnames(M))
  stopifnot(min(Q) >= 0L, max(Q) <= 255L)
  list(raw = matrix(as.raw(Q), nrow = nrow(Q), dimnames = dimnames(Q)), gmax = gmax)
}

# A long-format bin/gene/mean table becomes a dense [bin x gene] matrix. Bins
# absent for a gene are genuine zeros - the export drops them - so the matrix
# starts at zero rather than NA.
densify <- function(bin, gene, value, bins) {
  g <- sort(unique(gene))
  M <- matrix(0, nrow = length(bins), ncol = length(g), dimnames = list(NULL, g))
  i <- match(bin, bins); j <- match(gene, g)
  stopifnot(!anyNA(i), !anyNA(j))
  M[cbind(i, j)] <- value
  M
}

read_expr <- function(path, bins) {
  stopifnot(file.exists(path))
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    M <- readRDS(path)
    if (is.list(M) && !is.null(M$m)) M <- M$m
    stopifnot(is.matrix(M), is.numeric(M), !is.null(colnames(M)))
    stopifnot(nrow(M) == length(bins))
    say("  --expr: matrix %d bins x %s genes", nrow(M), format(ncol(M), big.mark = ","))
    M
  } else {
    x <- read.csv(path, stringsAsFactors = FALSE)
    need <- c("bin", "gene", "mean")
    if (!all(need %in% names(x)))
      stop("--expr csv needs columns bin, gene, mean; found: ", paste(names(x), collapse = ", "))
    say("  --expr: %s rows, %s genes", format(nrow(x), big.mark = ","),
        format(length(unique(x$gene)), big.mark = ","))
    densify(x$bin, x$gene, x$mean, bins)
  }
}

build_bundle <- function(DATA, CC = NULL, EXPR = NULL, CARRY = NULL,
                         GENES_OUT = NULL, CHUNK = 64L, allow_partial = FALSE) {
  stopifnot(dir.exists(DATA))
  rd <- function(...) {
    p <- file.path(DATA, ...)
    if (!file.exists(p)) stop("missing source file: ", p)
    read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  }

  # ---- 1. tables that ship as-is -----------------------------------------
  say("reading source tables from %s", DATA)
  B <- list(
    screen   = rd("Figure_5", "fig5A_top15.csv"),
    funnel   = rd("Figure_5", "fig5A_funnel.csv"),
    h2h      = rd("Figure_5", "fig5C_headtohead.csv"),
    ls_expr  = rd("Figure_5", "fig5D_ls.csv"),
    cyto     = rd("Figure_5", "fig5D_cyto.csv"),
    meta_ped = rd("Figure_5", "fig5F_meta_all.csv"),
    meta_ad  = rd("Figure_5", "fig5F_adult_beataml_meta.csv"),
    sc_ad    = rd("Figure_S23", "adult_scrna_overall.csv"),
    radar    = rd("Figure_6", "fig6D_radar.csv"),
    tox_hema = rd("Figure_6", "fig6tox_hema.csv"),
    tox_org  = rd("Figure_6", "fig6tox_organ.csv"),
    tox_ct   = (function(x) x[, c("gene", "label", "group", "cell_type",
                                  "pct", "mean", "n")])(
                 rd("Figure_6", "fig6tox_celltypes.csv")),
    combos   = rd("Figure_7", "fig7_allpairs.csv"),
    lsspec   = rd("Figure_S16", "funnel_ls_specific.csv"),
    lscomp   = rd("Figure_S16", "funnel_ls_specific_composition.csv"),
    ebins    = rd("Figure_2", "fig2A_bins.csv"))
  say("  screen %s genes | h2h %d | radar %d deep-profiled | ebins %s bins",
      format(nrow(B$screen), big.mark = ","), nrow(B$h2h), nrow(B$radar),
      format(nrow(B$ebins), big.mark = ","))
  stopifnot(nrow(B$ebins) > 0, !anyDuplicated(B$ebins$bin))

  # ---- 2. embedding, with the sample each cell came from ------------------
  emb <- rd("Figure_2", "fig2A_umap_cells.csv")
  key <- rd("Figure_4/v2", "tnk_subset_percell.csv")
  idx <- function(bc) sub("^.*-", "", bc)
  map <- unique(data.frame(i = idx(key$barcode), sample = key$SampleID,
                           stringsAsFactors = FALSE))
  if (anyDuplicated(map$i))
    stop("barcode suffix -> SampleID is not 1:1; cannot label samples")
  samp <- map$sample[match(idx(emb$barcode), map$i)]
  if (anyNA(samp))
    stop(sum(is.na(samp)), " cells have a barcode suffix absent from the ",
         "sample key; refusing to emit an embedding with unlabelled samples")
  LS_LEVELS <- paste0("LS_", sort(as.integer(sub("^LS_", "", unique(emb$LS[nzchar(emb$LS)])))))
  B$emb <- data.frame(
    UMAP1       = round(emb$UMAP1, 3),
    UMAP2       = round(emb$UMAP2, 3),
    compartment = factor(emb$compartment),
    cell_type   = factor(emb$cell_type),
    LS          = factor(emb$LS, levels = c("", LS_LEVELS)),
    bin         = as.integer(emb$bin),
    sample      = factor(samp),
    marrow      = factor(ifelse(grepl("^0_Healthy", samp),
                                "Healthy donor marrow", "AML marrow")))
  say("  emb %s cells | %d cell types | %d states | %d samples (%d donors)",
      format(nrow(B$emb), big.mark = ","), nlevels(B$emb$cell_type), length(LS_LEVELS),
      nlevels(B$emb$sample), sum(grepl("^0_Healthy", levels(B$emb$sample))))
  stopifnot(all(B$emb$bin %in% B$ebins$bin))

  # ---- 3. T/NK compartment, its CITE gates and its receptors --------------
  tnk <- rd("Figure_S12", "figS11_tnk_surface_percell.csv")
  mk  <- grep("_pos$", names(tnk), value = TRUE)
  if (!length(mk)) stop("no *_pos CITE-seq columns in figS11_tnk_surface_percell.csv")
  bn <- emb$bin[match(tnk$barcode, emb$barcode)]
  if (anyNA(bn))
    stop(sum(is.na(bn)), " T/NK cells have no barcode in the whole-atlas ",
         "embedding; cannot assign grid bins")
  SUBSET_LV <- names(sort(table(tnk$subset), decreasing = TRUE))
  B$tnk <- data.frame(UMAP1 = round(tnk$UMAP1, 3), UMAP2 = round(tnk$UMAP2, 3),
                      subset = factor(tnk$subset, levels = SUBSET_LV),
                      group = factor(tnk$group),
                      sample = factor(tnk$SampleID),
                      lapply(tnk[, mk, drop = FALSE], as.integer),
                      bin = as.integer(bn), check.names = FALSE)
  B$tnk_markers <- mk
  tc <- rd("Figure_4/v2", "tcr_percell_exh.csv")
  if (!all(c("subset", "group", "clone_size", "exhaustion", "expanded") %in% names(tc)))
    stop("tcr_percell_exh.csv is missing one of subset/group/clone_size/exhaustion/expanded")
  if (!all(unique(tc$expanded) %in% c("expanded", "unexpanded")))
    stop("unexpected values in tcr expanded: ",
         paste(setdiff(unique(tc$expanded), c("expanded", "unexpanded")), collapse = ", "))
  B$tcr <- data.frame(subset = factor(tc$subset, levels = SUBSET_LV),
                      group = factor(tc$group),
                      clone_size = as.integer(tc$clone_size),
                      exhaustion = round(as.numeric(tc$exhaustion), 4),
                      expanded = tc$expanded == "expanded")
  say("  tnk %s cells x %d markers | tcr %s cells, %s in an expanded clone",
      format(nrow(B$tnk), big.mark = ","), length(mk),
      format(nrow(B$tcr), big.mark = ","),
      format(sum(B$tcr$expanded), big.mark = ","))

  # ---- 4. differential tests ---------------------------------------------
  td <- rd("Figure_4/v2", "tnk_aml_vs_hbm_de_bySubset.csv")
  td$subset <- factor(td$subset, levels = SUBSET_LV)
  if (anyNA(td$subset)) stop("tnk_de carries a subset absent from the T/NK table")
  for (cl in c("logFC", "pct_AML", "pct_HBM")) td[[cl]] <- round(td[[cl]], 4)
  B$tnk_de <- td
  lsu <- rd("Figure_3/v3", "LS_upgenes_wilcoxon.csv")
  lsu$log2FC <- round(lsu$log2FC, 3)
  B$ls_up <- lsu
  say("  tnk_de %s rows, %s genes x %d subsets | ls_up %s rows, %s genes x %d states",
      format(nrow(B$tnk_de), big.mark = ","),
      format(length(unique(B$tnk_de$gene)), big.mark = ","),
      length(unique(B$tnk_de$subset)), format(nrow(B$ls_up), big.mark = ","),
      format(length(unique(B$ls_up$gene)), big.mark = ","),
      length(unique(B$ls_up$LS)))

  # ---- 5. bulk expression by cytogenetic subtype -------------------------
  cb <- rbind(
    transform(rd("Figure_S20", "figS_cyto_pediatric_violin.csv"),
              cohort = "TARGET (paediatric)"),
    transform(rd("Figure_S20", "figS_cyto_adult_violin.csv"),
              cohort = "Beat AML (adult)"))
  B$cyto_bulk <- data.frame(gene = cb$gene, cyto = cb$cyto,
                            log2cpm = round(cb$log2cpm, 4), cohort = cb$cohort,
                            stringsAsFactors = FALSE)
  ngb <- length(unique(B$cyto_bulk$gene))
  for (co in unique(B$cyto_bulk$cohort))
    say("  cyto_bulk %-20s %d patients x %d genes", co,
        sum(B$cyto_bulk$cohort == co) / ngb, ngb)

  # ---- 6. inlined expression for the screened genes ----------------------
  gx <- rd("Figure_2", "fig2A_gene_bin_expr.csv")
  M  <- densify(gx$bin, gx$gene, gx$mean, B$ebins$bin)
  say("  fig2A_gene_bin_expr %s rows -> %d bins x %s genes",
      format(nrow(gx), big.mark = ","), nrow(M), format(ncol(M), big.mark = ","))

  # ---- 7. the all-gene scope, and the runtime chunks ---------------------
  ALL <- NULL
  if (!is.null(EXPR)) {
    say("building the all-gene scope from --expr")
    ALL <- read_expr(EXPR, B$ebins$bin)
    # The six clinical comparators the shipped bundle inlined beyond the CSV.
    extra <- setdiff(c("ABCA7", "CD33", "CD38", "CLEC12A", "FLT3", "IL3RA"), colnames(M))
    extra <- intersect(extra, colnames(ALL))
    if (length(extra)) {
      M <- cbind(M, ALL[, extra, drop = FALSE])
      M <- M[, sort(colnames(M)), drop = FALSE]
      say("  inlined %d comparators alongside the screen genes: %s",
          length(extra), paste(extra, collapse = ", "))
    }
    q <- quantise(ALL)
    B$all_genes  <- colnames(ALL)
    B$all_gmax   <- q$gmax
    B$chunk_size <- as.integer(CHUNK)
    say("  all_genes %s | all_gmax %.4f | chunk_size %d",
        format(length(B$all_genes), big.mark = ","), B$all_gmax, B$chunk_size)
    if (!is.null(GENES_OUT)) {
      dir.create(GENES_OUT, showWarnings = FALSE, recursive = TRUE)
      nch <- ceiling(ncol(q$raw) / CHUNK)
      for (ci in seq_len(nch) - 1L) {
        j <- (ci * CHUNK + 1L):min((ci + 1L) * CHUNK, ncol(q$raw))
        saveRDS(q$raw[, j, drop = FALSE],
                file.path(GENES_OUT, sprintf("g%03d.rds", ci)), compress = "gzip")
      }
      writeLines(paste0("[", paste0('"', B$all_genes, '"', collapse = ","), "]"),
                 file.path(GENES_OUT, "index.json"))
      say("  wrote %d chunks + index.json to %s", nch, GENES_OUT)
    } else {
      say("  NOTE: --genes-out not given, so the runtime chunks were NOT written.")
    }
  } else if (!is.null(CARRY)) {
    stopifnot(file.exists(CARRY))
    old <- readRDS(CARRY)
    need <- c("all_genes", "all_gmax", "chunk_size")
    if (!all(need %in% names(old))) stop("--carry bundle lacks ", paste(need, collapse = "/"))
    B[need] <- old[need]
    CARRIED_GMAT <- old
    say("CARRIED OVER from %s, not rebuilt: all_genes (%s), all_gmax (%.4f), chunk_size (%d).",
        CARRY, format(length(B$all_genes), big.mark = ","), B$all_gmax, B$chunk_size)
    say("  The genes/ chunk files on disk are assumed to match and were not touched.")
  } else {
    msg <- paste("no --expr and no --carry: the all-gene scope would be absent and",
                 "the atlas would offer only the screened genes")
    if (!allow_partial) stop(msg, ". Pass --allow-partial to build anyway.")
    say("WARNING: %s.", msg)
  }
  qm <- quantise(M)
  B$gmat <- qm$raw; B$gmax <- qm$gmax; B$genes_mapped <- colnames(M)
  # The CSV holds the screen genes only. When --expr is unavailable but --carry
  # is given, bring the extra inlined columns across rather than letting them
  # fall back to a runtime fetch. Only safe because both are quantised against
  # the same maximum, which is asserted rather than assumed.
  if (exists("CARRIED_GMAT", inherits = FALSE)) {
    extra <- setdiff(colnames(CARRIED_GMAT$gmat), colnames(B$gmat))
    if (length(extra)) {
      if (!isTRUE(all.equal(CARRIED_GMAT$gmax, B$gmax, tolerance = 1e-12)))
        stop("cannot carry inlined columns: the carried bundle quantised gmat ",
             "against ", CARRIED_GMAT$gmax, " and this build against ", B$gmax)
      shared <- intersect(colnames(CARRIED_GMAT$gmat), colnames(B$gmat))
      if (!identical(CARRIED_GMAT$gmat[, shared, drop = FALSE],
                     B$gmat[, shared, drop = FALSE]))
        stop("cannot carry inlined columns: the shared genes do not match the ",
             "carried bundle, so the two were not built from the same export")
      B$gmat <- cbind(B$gmat, CARRIED_GMAT$gmat[, extra, drop = FALSE])
      B$gmat <- B$gmat[, sort(colnames(B$gmat)), drop = FALSE]
      B$genes_mapped <- colnames(B$gmat)
      say("  carried %d inlined columns absent from the CSV: %s",
          length(extra), paste(extra, collapse = ", "))
    }
  }
  say("  gmat %d bins x %s genes | gmax %.4f", nrow(B$gmat),
      format(ncol(B$gmat), big.mark = ","), B$gmax)

  # ---- 8. cell-cell communication ----------------------------------------
  if (!is.null(CC)) {
    stopifnot(dir.exists(CC))
    norm <- sort(setdiff(levels(B$emb$cell_type), LEUK_TYPES))
    B$cc_leuk   <- c(LS_LEVELS, "AML (unassigned)")
    B$cc_normal <- norm
    B$cc_groups <- c(B$cc_leuk, B$cc_normal)
    say("building cell-cell communication over %d groups (%d leukemic, %d normal)",
        length(B$cc_groups), length(B$cc_leuk), length(B$cc_normal))
    files <- file.path(CC, sprintf("g%02d.rds", seq_along(B$cc_groups) - 1L))
    miss <- files[!file.exists(files)]
    if (length(miss))
      stop(length(miss), " of ", length(files), " CellChat group files missing, ",
           "first: ", miss[1], ". Refusing to build cc_pairs from a partial run.")
    pool <- do.call(rbind, lapply(files, function(f)
      readRDS(f)[, c("source", "target", "name", "prob")]))
    # Each group file holds every link touching that group, as source or as
    # target, so the pool carries all but the self-links twice.
    id <- paste(pool$source, pool$target, pool$name, sep = "\r")
    d  <- pool[!duplicated(id), ]
    say("  %s pooled rows -> %s distinct significant links",
        format(nrow(pool), big.mark = ","), format(nrow(d), big.mark = ","))
    k <- paste(d$source, d$target, sep = "\r")
    s <- tapply(d$prob, k, sum); n <- tapply(d$prob, k, length)
    pt <- do.call(rbind, strsplit(names(s), "\r", fixed = TRUE))
    cp <- data.frame(source = factor(pt[, 1], levels = B$cc_groups),
                     target = factor(pt[, 2], levels = B$cc_groups),
                     strength = as.numeric(s), n = as.integer(n))
    if (anyNA(cp$source) || anyNA(cp$target))
      stop("a CellChat link names a group outside the ", length(B$cc_groups),
           " this atlas defines")
    B$cc_pairs <- cp[order(as.integer(cp$source), as.integer(cp$target)), ]
    rownames(B$cc_pairs) <- NULL
    say("  cc_pairs %s ordered pairs", format(nrow(B$cc_pairs), big.mark = ","))

    # interaction -> gene, complexes expanded. From CellChatDB when it is
    # installed, otherwise from a cc_gene_map.csv deposited beside the run.
    side <- file.path(CC, "cc_gene_map.csv")
    if (requireNamespace("CellChat", quietly = TRUE)) {
      db  <- CellChat::CellChatDB.human
      cx  <- db$complex
      expand <- function(u) {
        u <- u[nzchar(u)]
        out <- unlist(lapply(u, function(x)
          if (x %in% rownames(cx)) unlist(cx[x, ], use.names = FALSE) else x))
        unique(out[nzchar(out) & !is.na(out)])
      }
      ia <- db$interaction
      gm <- do.call(rbind, lapply(seq_len(nrow(ia)), function(i) {
        g <- expand(c(ia$ligand[i], ia$receptor[i]))
        if (!length(g)) return(NULL)
        data.frame(interaction_name = ia$interaction_name[i], gene = g,
                   stringsAsFactors = FALSE)
      }))
      B$cc_gene_map <- unique(gm)   # see the note above about the 89 duplicates
      say("  cc_gene_map from CellChatDB.human: %s rows, %s interactions, %s genes",
          format(nrow(B$cc_gene_map), big.mark = ","),
          format(length(unique(B$cc_gene_map$interaction_name)), big.mark = ","),
          format(length(unique(B$cc_gene_map$gene)), big.mark = ","))
      write.csv(B$cc_gene_map, side, row.names = FALSE)
      say("  wrote %s so a later run needs no CellChat install", side)
    } else if (file.exists(side)) {
      B$cc_gene_map <- read.csv(side, stringsAsFactors = FALSE)
      say("  cc_gene_map from %s (CellChat not installed): %s rows", side,
          format(nrow(B$cc_gene_map), big.mark = ","))
    } else {
      stop("cc_gene_map needs either the CellChat package or ", side,
           ". It cannot be recovered from the link files: the map covers every ",
           "database interaction, not only the observed ones, and it expands ",
           "complexes that the interaction name does not spell out.")
    }
    B$cc_db_genes <- sort(unique(B$cc_gene_map$gene))
  } else {
    msg <- "no --cc: the Cell talk tab would have no data"
    if (!allow_partial) stop(msg, ". Pass --allow-partial to build anyway.")
    say("WARNING: %s.", msg)
  }

  # ---- 9. restrict de_all to the atlas's own genes ------------------------
  de <- rd("Figure_2/v2_recluster", "group_DE_wholecohort.csv")
  keep <- if (!is.null(B$all_genes)) B$all_genes else B$genes_mapped
  de <- de[de$gene %in% keep, ]
  de$log2FC <- round(de$log2FC, 4)
  B$de_all <- de
  say("  de_all %s of %s tested genes kept (those the atlas quantifies)",
      format(nrow(B$de_all), big.mark = ","),
      format(nrow(rd("Figure_2/v2_recluster", "group_DE_wholecohort.csv")), big.mark = ","))

  # app.R colours these categories by NAME, through the KLASS vector. A category
  # renamed in the source data simply misses the palette and draws uncoloured,
  # with no error, so print the labels on every build and check them against
  # KLASS. This has already happened once: the bundle that shipped on
  # 2026-10-21 said "ppac" for CD7, TNFRSF4, ABCA7 and ITGAX where these CSVs
  # now say "state_specific", so KLASS still has a ppac key. A rebuild from the
  # current source needs that key renamed in the same commit.
  cats <- sort(unique(c(B$radar$group, B$h2h$klass)))
  say("\ncategory labels in this build (app.R must colour each by name):\n  %s",
      paste(cats, collapse = ", "))

  B
}

# Diff two bundles object by object. Used by --check; this is what proves a
# rebuild reproduces what shipped, rather than merely running without error.
compare_bundles <- function(new, ref, tol = 1e-9) {
  nm <- union(names(ref), names(new)); bad <- 0L
  for (n in nm) {
    if (!n %in% names(new)) { say("  %-14s MISSING from the rebuild", n); bad <- bad + 1L; next }
    if (!n %in% names(ref)) { say("  %-14s extra in the rebuild", n); next }
    a <- new[[n]]; b <- ref[[n]]
    chr <- function(d) lapply(d, function(z) if (is.factor(z)) as.character(z) else z)
    same <- if (is.data.frame(a) && is.data.frame(b)) {
      identical(dim(a), dim(b)) && identical(names(a), names(b)) &&
        isTRUE(all.equal(chr(a), chr(b), tolerance = tol))
    } else isTRUE(all.equal(a, b, tolerance = tol))
    if (same) { say("  %-14s OK", n); next }
    # Same rows in a different order is not a data difference. Say which it is,
    # so a cosmetic reordering cannot be mistaken for a changed value.
    if (is.data.frame(a) && is.data.frame(b) && identical(dim(a), dim(b)) &&
        identical(names(a), names(b))) {
      ka <- do.call(paste, c(chr(a), sep = "\r"))
      kb <- do.call(paste, c(chr(b), sep = "\r"))
      if (setequal(ka, kb)) { say("  %-14s same rows, different order", n); next }
    }
    say("  %-14s *** DIFFERS", n); bad <- bad + 1L
  }
  say("%d of %d objects differ", bad, length(nm))
  invisible(bad)
}

if (!interactive() && sys.nframe() == 0L) {
  a <- commandArgs(TRUE)
  opt <- function(f, d = NULL) { i <- match(f, a); if (is.na(i)) d else a[i + 1L] }
  pos <- a[!grepl("^--", a)]
  pos <- pos[!pos %in% unlist(lapply(c("--cc","--expr","--genes-out","--carry","--out","--check"),
                                     function(f) opt(f)))]
  src <- if (length(pos)) pos[1] else "source_data"
  out <- opt("--out", file.path("app", "bundle.rds"))
  B <- build_bundle(src, CC = opt("--cc"), EXPR = opt("--expr"),
                    CARRY = opt("--carry"), GENES_OUT = opt("--genes-out"),
                    allow_partial = "--allow-partial" %in% a)
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  # gzip, not xz: xz saves ~0.5 MB and decompresses about 5x slower, and
  # startup time is what the browser feels.
  saveRDS(B, out, compress = "gzip")
  say("\nbundle -> %s (%.1f MB), %d objects", out, file.size(out) / 1e6, length(B))
  chk <- opt("--check")
  if (!is.null(chk)) { say("\nchecking against %s", chk); compare_bundles(B, readRDS(chk)) }
}
