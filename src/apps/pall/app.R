# Pediatric ALL treatment resistance browser
# ---------------------------------------------------------------------------
# Pang, Prieto et al.: "Single-Cell Sequencing Reveals Extensive Genetic
# Diversity Underlying Pediatric ALL Treatment Complexity". Every number here
# comes from the manuscript's figure-level source data; nothing is recomputed
# from raw sequence and nothing is simulated.
#
#   Rscript src/apps/pall/prep_bundle.R SOURCE_DATA bundle.rds
#   shiny::runApp("src/apps/pall")
# ---------------------------------------------------------------------------

library(shiny); library(bslib); library(ggplot2); library(DT); library(ape); library(scales)

B <- local({
  for (p in c("bundle.rds", file.path("..", "bundle.rds"),
              file.path("src", "apps", "pall", "bundle.rds")))
    if (file.exists(p)) return(readRDS(p))
  stop("bundle.rds not found - build it with prep_bundle.R first")
})
burden <- B$burden; pan <- B$pan; ras <- B$ras
drug <- B$drug; DMAT <- B$dmat; sj <- B$sj
cells <- B$cells; GMAT <- B$gmat; EMERGENT <- B$emergent
rec <- B$rec; am <- B$am
TREES <- B$trees; TIPS <- B$tips; AF <- B$af
CELLMETA <- B$cellmeta; SIG <- B$sig

# Newick is parsed here rather than shipped as a serialised tree, so the bundle
# stays in kilobytes. Labels come through as 4295_A10 or 4295.F1 depending on the
# file; the annotation table is keyed on the normalised 4295-A10 form.
TREE_CHOICES <- local({
  n <- vapply(TREES, function(x) length(gregexpr("[,(]", x)[[1]]) , integer(1))
  lab <- c("4295" = "Patient 4295", "445" = "Patient 445", "417" = "Patient 417",
           "4084" = "Patient 4084", "Invitro" = "In vitro benchmark",
           "clone" = "Clone tree, both timepoints")
  k <- intersect(names(lab), names(TREES))
  setNames(k, sprintf("%s", lab[k]))
})
norm_tip <- function(x) gsub("[._]", "-", x)

CARD <- "#8C1515"; TEAL <- "#2F5D70"; AMBER <- "#F6A30C"; GREY <- "#B9C0C7"
GRP <- c(MRD_candidate = "#D55E00", control_relapse = "#0072B2", control_founding = "#E69F00")
GRP_LAB <- c(MRD_candidate = "MRD candidate", control_relapse = "Resistance-associated",
             control_founding = "Known cancer initiator")
COND_COL <- c(Bulk = GREY, DMSO = TEAL, `Pred-Hi` = CARD, `DNR-Hi` = AMBER)
MARKERS <- c("CD19_M", "CD34_M", "CD10_M", "CD20_M")

theme_lab <- function(rot = 0) {
  theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          axis.text.x = element_text(angle = rot, hjust = if (rot > 0) 1 else .5))
}
note <- function(...) div(class = "text-muted", style = "font-size:.85rem;margin-top:.5rem", ...)
ext  <- function(l, u) a(l, href = u, target = "_blank", rel = "noopener")

DRUG_GENES <- sort(unique(sub(" .*$", "", colnames(DMAT))))
ALL_GENES  <- sort(unique(c(DRUG_GENES, rec$gene, am$gene, sj$gene)))

# --- ui --------------------------------------------------------------------
ui <- page_navbar(
  id = "nav",
  title = "Pediatric ALL Treatment Resistance",
  # Panels scroll instead of filling. With fill on, bslib divides the viewport
  # height between the cards in a panel, so a tab with five cards rendered every
  # plot 17 pixels tall no matter what height plotOutput asked for.
  fillable = FALSE,
  theme = bs_theme(version = 5, primary = "#2F5D70"),
  header = tags$style(HTML(
    ".navbar .navbar-brand{font-size:1.45rem;font-weight:700}",
    ".navbar .nav-link{font-size:1.18rem;font-weight:600;padding:.5rem 1rem}",
    ".navbar .nav-link.active{font-weight:700}")),
  sidebar = sidebar(
    width = 300,
    conditionalPanel("input.nav == 'Drug response'",
      selectInput("dpat", "Patient", sort(unique(drug$patient))),
      selectizeInput("dmut", "Highlight a mutation", choices = colnames(DMAT),
                     options = list(maxOptions = 200))),
    conditionalPanel("input.nav == 'Single cells'",
      selectInput("cpat", "Patient", sort(unique(CELLMETA$patient))),
      radioButtons("cfill", "Colour the paired-sample plot by",
                   c("Clone" = "clone", "Timepoint" = "timepoint",
                     "Chromosome 4 deletion" = "Chr4_Deletion",
                     "Chromosome 6 deletion" = "Chr6_1_Deletion"))),
    conditionalPanel("input.nav == 'Phylogeny'",
      selectInput("tree", "Tree", TREE_CHOICES),
      radioButtons("ttype", "Layout",
                   c("Phylogram" = "phylogram", "Fan" = "fan", "Unrooted" = "unrooted")),
      uiOutput("tipcol_ui"),
      sliderInput("bootmin", "Label bootstrap support of at least", 0, 100, 50, step = 5),
      checkboxInput("tiplab", "Show cell labels", FALSE)),
    conditionalPanel("input.nav == 'Genes'",
      selectizeInput("gene", "Highlight a gene", choices = sort(unique(rec$gene)),
                     selected = if ("TBL1XR1" %in% rec$gene) "TBL1XR1" else sort(unique(rec$gene))[1],
                     options = list(maxOptions = 200))),
    hr(),
    note(strong("Pang, Prieto ", em("et al."), "."),
         " Single-cell sequencing reveals extensive genetic diversity underlying ",
         "pediatric ALL treatment complexity. Every panel is drawn from the ",
         "manuscript's figure-level source data.")
  ),

  nav_panel("Hidden diversity",
    layout_columns(col_widths = c(6, 6),
      card(card_header("Bulk sequencing misses most of the mutations"),
           plotOutput("burden_plot", height = 400), uiOutput("burden_head"),
           note("Whole-genome calls from the five patients sequenced both ways, each ",
                "corrected for that sample's own detection sensitivity. A mutation private ",
                "to one clone is diluted below the detection floor of a bulk sample, which ",
                "is why the single cells sit so far above it. These are the Figure 3 ",
                "whole-genome numbers; the exome experiment in Figure 2 measures the same ",
                "effect on a different scale.")),
      card(card_header("Where pediatric ALL sits among childhood cancers"),
           plotOutput("pan_plot", height = 400),
           note("SNVs per megabase in a published survey of paediatric tumours. ",
                "ALL sits at the bottom of the range, which is the observation this ",
                "study set out to reconcile with the complexity of its treatment."))),
    card(card_header("Per-patient detail"), DTOutput("burden_tbl"))),

  nav_panel("RAS",
    card(card_header("Activating RAS mutations found by error-corrected sequencing"),
         plotOutput("ras_plot", height = 420),
         note("Each point is one activating mutation. Bulk sequencing reported a single ",
              "RAS mutation in each of these patients; error-corrected sequencing finds ",
              "several more at low allele frequency, which is why they were missed.")),
    layout_columns(col_widths = c(7, 5),
      card(card_header("Allele frequency by codon"), plotOutput("ras_codon", height = 340)),
      card(card_header("Every RAS mutation"), DTOutput("ras_tbl")))),

  nav_panel("Drug response",
    card(card_header(textOutput("drug_title")), plotOutput("drug_heat", height = 460),
         note("Mutant allele frequency in percent, one row per sequenced sample and one ",
              "column per recurrent mutation. Pred-Hi is prednisolone and DNR-Hi is ",
              "daunorubicin, each against its own DMSO control and the diagnostic bulk ",
              "sample. A column that rises under one drug and not the other marks a ",
              "population with differential sensitivity.")),
    layout_columns(col_widths = c(6, 6),
      card(card_header(textOutput("dmut_title")), plotOutput("dmut_plot", height = 340)),
      card(card_header("SJETV077 across nine ex vivo conditions"),
           plotOutput("sj_plot", height = 340),
           note("A separate single-patient experiment covering six agents plus controls.")))),

  nav_panel("Single cells",
    card(card_header(textOutput("qc_title")),
      layout_columns(col_widths = c(6, 6),
        plotOutput("qc_plot", height = 380), plotOutput("sig_plot", height = 380)),
      uiOutput("qc_head"),
      note("Every patient's cells are here. Allelic dropout and depth are the two ",
           "measurements that decide whether a single-cell genome can be called at all, ",
           "and the signature panel is the COSMIC exposure fitted to each cell's own ",
           "mutations. The panels below need the paired pretreatment and post-treatment ",
           "samples, which exist for one patient only.")),

    h4(textOutput("paired_head"), style = "margin:1.6rem 0 .6rem;font-weight:700"),

    layout_columns(col_widths = c(7, 5),
      card(card_header("115 single-cell genomes, before and after treatment"),
           plotOutput("cell_plot", height = 430),
           note("Two paired samples from the same patient: 4272 before treatment and ",
                "4295 after. Axes are the measured surface-marker intensities used to ",
                "separate leukemic from normal and premalignant cells.")),
      card(card_header("Clone composition"), plotOutput("clone_plot", height = 430),
           note("Each bar is the percentage of that timepoint's cells, not a raw count, ",
                "because 30 cells were sequenced before treatment against 85 after. Clones ",
                strong(paste(EMERGENT, collapse = " and ")),
                " are absent from the pretreatment sample entirely and appear only after ",
                "treatment, but they are 3 and 1 cells, so read them as the observation ",
                "they are rather than as a reliable frequency."))),
    card(card_header("Allele frequency before and after treatment"),
      plotOutput("af_plot", height = 470),
      uiOutput("af_head"),
      note("Each point is one somatic variant, pooling the alt and total reads of every cell ",
           "in that sample, so this is a pseudobulk allele frequency rather than a per-cell ",
           "call. Points above the diagonal rose under treatment. Read counts come from a ",
           "targeted panel of 31 variants, so this is not a genome-wide survey, and the ",
           "pretreatment sample carries 30 cells against 85 after, which makes the ",
           "pretreatment estimate the noisier of the two.")),

    card(card_header("Every variant, before and after"), DTOutput("af_tbl")),

    card(card_header("Genotypes across 31 variants"), plotOutput("geno_plot", height = 420),
         note("Presence or absence of each somatic variant in each cell, cells ordered by clone.")),

    card(card_header(textOutput("pie_title")),
      layout_columns(col_widths = c(5, 7),
        plotOutput("sig_pie", height = 430), plotOutput("sig_bar", height = 430)),
      uiOutput("pie_head"),
      note("This pie pools every mutation from every cell of the patient and splits that total ",
           "by signature. It is not an average of the per-cell percentages: a cell that ",
           "contributed 12,000 mutations counts for more than one that contributed 300, which is ",
           "what treating each cell as one unit would wrongly do. The per-cell pies on the ",
           "Phylogeny tab are the other view, each one that single cell's own mutations split by ",
           "signature. Signatures under 3% are pooled here and under 2% in the bars. The in vitro ",
           "benchmark is a cell line rather than a patient, and is left out of the four-patient ",
           "comparison."))),

  nav_panel("Phylogeny",
    card(fill = FALSE, card_header(textOutput("tree_title")),
         plotOutput("tree_plot", height = "auto"),
         uiOutput("tree_legend"), uiOutput("tree_head"),
         note("A tip pie is that one cell's own mutations split by signature. Where a tree ",
              "holds both samples, the branch leading to each cell is coloured by the sample it ",
              "came from, so before and after treatment read off the topology directly. ",
              "The patient trees are maximum-likelihood phylogenies built by CellPhy from ",
              "somatic single-nucleotide variants, with support from 100 bootstrap replicates; ",
              "branch lengths are substitutions per site. The clone tree is the ConDoR ",
              "topology behind Figure 7A, which carries no branch lengths, so it is drawn as ",
              "a cladogram and the fan and unrooted layouts show topology only."))),

  nav_panel("Genes",
    card(card_header("Recurrence against predicted pathogenicity"),
         plotOutput("gene_scatter", height = 440),
         uiOutput("gene_count"),
         note("Horizontal axis is recurrence above what the gene's coding length predicts, on a ",
              "log scale; vertical axis is the mean AlphaMissense score of the missense variants ",
              "observed in it. Upper right is a gene that is both hit more often than ",
              "expected and predicted damaging. The two axes come from different cohort ",
              "sets, so a gene's pair of values is not derived from identical samples.")),
    layout_columns(col_widths = c(7, 5),
      card(card_header(textOutput("am_title")), plotOutput("am_plot", height = 360),
           note("Dashed lines are the published AlphaMissense thresholds: likely benign ",
                "below 0.34, likely pathogenic above 0.564.")),
      card(card_header("Gene table"), DTOutput("gene_tbl"))))
)

# --- server ----------------------------------------------------------------
server <- function(input, output, session) {

  # ---- Hidden diversity ----
  # There is one bulk sample but several single cells per patient, so the cells
  # are shown as a distribution rather than joined to the bulk point by a line.
  output$burden_plot <- renderPlot({
    d <- burden; d$assay <- factor(d$assay, levels = c("Bulk", "Single"))
    ggplot(d, aes(patient, corrected_som_per_mb, colour = assay)) +
      geom_point(position = position_jitterdodge(jitter.width = .25, dodge.width = .7, seed = 5),
                 size = 3, alpha = .9) +
      stat_summary(fun = median, geom = "crossbar", width = .45,
                   position = position_dodge(width = .7), linewidth = .4, show.legend = FALSE) +
      scale_colour_manual(values = c(Bulk = "#4A555F", Single = CARD), name = NULL) +
      labs(x = NULL, y = "sensitivity-corrected somatic mutations per Mb") + theme_lab()
  })

  # stated from the table in front of the reader, not quoted from the paper
  output$burden_head <- renderUI({
    b <- burden$corrected_som_per_mb[burden$assay == "Bulk"]
    s <- burden$corrected_som_per_mb[burden$assay == "Single"]
    note(sprintf(paste("Across these five patients the median single cell carries %.1f times the",
                       "corrected mutation density of its matched bulk sample (%.2f against %.2f",
                       "per Mb, %d cells and %d bulk samples)."),
                 median(s) / median(b), median(s), median(b), length(s), length(b)))
  })

  output$pan_plot <- renderPlot({
    d <- pan
    med <- stats::aggregate(snv_per_mb ~ cancer_type, d, median)
    d$cancer_type <- factor(d$cancer_type, levels = med$cancer_type[order(med$snv_per_mb)])
    d$is_all <- grepl("ALL", d$cancer_type)
    ggplot(d, aes(snv_per_mb, cancer_type, colour = is_all)) +
      geom_point(position = position_jitter(height = .18, seed = 1), size = 1.1, alpha = .6) +
      scale_x_continuous(trans = "log10") +
      scale_colour_manual(values = c(`FALSE` = GREY, `TRUE` = CARD), guide = "none") +
      labs(x = "SNVs per Mb (log scale)", y = NULL) + theme_lab()
  })

  output$burden_tbl <- renderDT({
    d <- burden[, c("patient", "assay", "total_somatic", "unique", "shared",
                    "sensitivity", "corrected_som_per_mb", "corrected_total")]
    names(d) <- c("Patient", "Assay", "Total somatic", "Unique", "Shared",
                  "Sensitivity", "Corrected som./Mb", "Corrected total")
    datatable(d, rownames = FALSE, options = list(pageLength = 10, dom = "tip")) |>
      formatRound(c("Sensitivity", "Corrected som./Mb"), 2) |>
      formatRound("Corrected total", 0)
  })

  # ---- RAS ----
  output$ras_plot <- renderPlot({
    d <- ras
    n <- stats::aggregate(AF ~ Patient, d, length); names(n)[2] <- "n"
    d$Patient <- factor(d$Patient, levels = n$Patient[order(-n$n)])
    ggplot(d, aes(AF, Patient, colour = Ras)) +
      geom_point(size = 3, alpha = .85) +
      scale_x_continuous(trans = "log10", labels = function(x) paste0(x * 100, "%")) +
      scale_colour_manual(values = c(KRAS = CARD, NRAS = TEAL), name = NULL) +
      labs(x = "mutant allele frequency (log scale)", y = NULL) + theme_lab()
  })

  output$ras_codon <- renderPlot({
    d <- ras; d$Location <- factor(d$Location, levels = sort(unique(d$Location)))
    ggplot(d, aes(Location, AF, colour = Ras)) +
      geom_point(position = position_jitter(width = .12, seed = 2), size = 3, alpha = .85) +
      scale_y_continuous(trans = "log10", labels = function(x) paste0(x * 100, "%")) +
      scale_colour_manual(values = c(KRAS = CARD, NRAS = TEAL), name = NULL) +
      labs(x = "codon", y = "allele frequency") + theme_lab()
  })

  output$ras_tbl <- renderDT({
    d <- ras[order(-ras$AF), c("Patient", "Ras", "Location", "AA_Change", "AF")]
    names(d) <- c("Patient", "Gene", "Codon", "Change", "Allele frequency")
    datatable(d, rownames = FALSE, options = list(pageLength = 8, dom = "ftip")) |>
      formatPercentage("Allele frequency", 2)
  })

  # ---- Drug response ----
  dsub <- reactive({ i <- drug$patient == input$dpat; list(meta = drug[i, ], m = DMAT[i, , drop = FALSE]) })

  output$drug_title <- renderText(sprintf("Patient %s: ex vivo response of every recurrent mutation", input$dpat))
  output$drug_heat <- renderPlot({
    s <- dsub(); m <- s$m
    keep <- colSums(m > 0, na.rm = TRUE) > 0
    validate(need(any(keep), "No mutations recorded for this patient."))
    m <- m[, keep, drop = FALSE]
    ord <- order(colMeans(m, na.rm = TRUE), decreasing = TRUE)
    m <- m[, ord, drop = FALSE]
    lab <- paste(s$meta$condition, s$meta$replicate)
    d <- data.frame(sample = factor(rep(lab, ncol(m)), levels = rev(lab[order(s$meta$condition)])),
                    mutation = factor(rep(colnames(m), each = nrow(m)), levels = colnames(m)),
                    af = as.vector(m), stringsAsFactors = FALSE)
    ggplot(d, aes(mutation, sample, fill = af)) +
      geom_tile(colour = "white", linewidth = .25) +
      scale_fill_viridis_c(option = "rocket", direction = -1, name = "AF %") +
      labs(x = NULL, y = NULL) + theme_lab(90) +
      theme(axis.text.x = element_text(size = 7), panel.grid = element_blank())
  })

  output$dmut_title <- renderText(sprintf("%s across conditions", input$dmut))
  output$dmut_plot <- renderPlot({
    req(input$dmut %in% colnames(DMAT))
    d <- drug; d$af <- DMAT[, input$dmut]
    d$condition <- factor(d$condition, levels = names(COND_COL))
    ggplot(d, aes(condition, af, colour = condition)) +
      geom_point(position = position_jitter(width = .12, seed = 3), size = 3, alpha = .9) +
      facet_wrap(~ patient, nrow = 1) +
      scale_colour_manual(values = COND_COL, guide = "none") +
      labs(x = NULL, y = "mutant allele frequency (%)") + theme_lab(45)
  })

  output$sj_plot <- renderPlot({
    d <- sj
    top <- names(sort(tapply(d$af, d$Mutation, max, na.rm = TRUE), decreasing = TRUE))[1:25]
    d <- d[d$Mutation %in% top, ]
    d$Mutation <- factor(d$Mutation, levels = rev(top))
    ggplot(d, aes(Treatment, Mutation, fill = af)) +
      geom_tile(colour = "white", linewidth = .25) +
      scale_fill_viridis_c(option = "mako", direction = -1, name = "AF %") +
      scale_y_discrete(labels = function(x) sub(":.*$", "", x)) +
      labs(x = NULL, y = NULL) + theme_lab(45) +
      theme(axis.text.y = element_text(size = 7), panel.grid = element_blank())
  })

  # ---- Single cells ----
  output$cell_plot <- renderPlot({
    d <- cells
    d$fill <- switch(input$cfill,
      clone = ifelse(is.na(d$clone), "unassigned", d$clone),
      timepoint = as.character(d$timepoint),
      ifelse(is.na(d[[input$cfill]]), "unknown",
             ifelse(d[[input$cfill]] == 1, "present", "absent")))
    d <- d[is.finite(d$CD19_M) & is.finite(d$CD34_M), ]
    validate(need(nrow(d) > 0, "No cells with both markers measured."))
    ggplot(d, aes(CD19_M, CD34_M, colour = fill, shape = timepoint)) +
      geom_point(size = 3, alpha = .85) +
      scale_colour_manual(values = grDevices::hcl.colors(length(unique(d$fill)), "Spectral"),
                          name = NULL) +
      scale_shape_manual(values = c(Pretreatment = 1, `Post-treatment` = 16), name = NULL) +
      labs(x = "CD19", y = "CD34") + theme_lab()
  })

  # Raw counts would mislead: 30 cells were sequenced before treatment against 85
  # after, so every clone looks larger afterwards. Each bar is the share of the
  # cells sequenced at that timepoint.
  output$clone_plot <- renderPlot({
    d <- cells; d$clone <- ifelse(is.na(d$clone), "unassigned", d$clone)
    tab <- as.data.frame(table(clone = d$clone, timepoint = d$timepoint))
    tot <- tapply(tab$Freq, tab$timepoint, sum)
    tab$pct <- tab$Freq / tot[as.character(tab$timepoint)] * 100
    tab$emergent <- tab$clone %in% EMERGENT
    lv <- levels(droplevels(tab$timepoint))
    labs_tp <- setNames(sprintf("%s (n = %d)", lv, tot[lv]), lv)
    tab$timepoint <- factor(labs_tp[as.character(tab$timepoint)], levels = labs_tp)
    ggplot(tab, aes(pct, clone, fill = timepoint)) +
      geom_col(position = "dodge") +
      geom_text(data = tab[tab$emergent & tab$pct > 0, ],
                aes(label = "emergent"), hjust = -.15, size = 3, colour = CARD) +
      scale_fill_manual(values = setNames(c(GREY, CARD), labs_tp), name = NULL) +
      scale_x_continuous(expand = expansion(mult = c(0, .28)),
                         labels = function(x) paste0(x, "%")) +
      labs(x = "share of the cells sequenced at that timepoint", y = NULL) + theme_lab()
  })

  output$geno_plot <- renderPlot({
    ord <- order(cells$timepoint, ifelse(is.na(cells$clone), "zz", cells$clone))
    ids <- cells$Index[ord]
    ids <- ids[ids %in% rownames(GMAT)]
    m <- GMAT[ids, , drop = FALSE]
    d <- data.frame(cell = factor(rep(ids, ncol(m)), levels = ids),
                    variant = factor(rep(colnames(m), each = nrow(m)), levels = colnames(m)),
                    present = as.vector(m) == 1)
    ggplot(d, aes(cell, variant, fill = present)) +
      geom_tile() +
      scale_fill_manual(values = c(`TRUE` = TEAL, `FALSE` = "#EDF0F2"), guide = "none") +
      scale_y_discrete(labels = function(x) sub("^.*_", "", x)) +
      labs(x = "cells, ordered by timepoint then clone", y = NULL) +
      theme_lab() + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 7),
                          panel.grid = element_blank())
  })

  # ---- per-patient cell quality and signatures ----
  cm <- reactive(CELLMETA[CELLMETA$patient == input$cpat, ])

  output$qc_title <- renderText(sprintf("Patient %s: %d single-cell genomes",
                                        input$cpat, nrow(cm())))
  output$paired_head <- renderText(
    sprintf("Paired pretreatment and post-treatment samples \u2014 patient 4295 only%s",
            if (identical(input$cpat, "4295")) "" else
              sprintf(" (you have patient %s selected above)", input$cpat)))

  output$qc_plot <- renderPlot({
    d <- cm(); d <- d[is.finite(d$ado) & is.finite(d$depth), ]
    validate(need(nrow(d) > 0, sprintf(
      "No allelic-dropout or depth measurements were reported for patient %s.", input$cpat)))
    ggplot(d, aes(depth, ado)) +
      geom_point(aes(size = n_snv), colour = CARD, alpha = .8) +
      scale_size_continuous(range = c(2, 7), name = "SNVs called", labels = scales::comma) +
      labs(x = "sequencing depth (fold)", y = "allelic dropout (%)") + theme_lab()
  })

  output$sig_plot <- renderPlot({
    d <- cm(); k <- intersect(d$cell, rownames(SIG))
    validate(need(length(k) > 0, sprintf("No signature fit was reported for patient %s.", input$cpat)))
    m <- SIG[k, , drop = FALSE]
    keep <- colSums(m) > 0
    tot <- sort(colSums(m[, keep, drop = FALSE]), decreasing = TRUE)
    top <- names(tot)[seq_len(min(8, length(tot)))]
    dd <- data.frame(signature = factor(rep(top, each = nrow(m)), levels = rev(top)),
                     exposure = as.vector(m[, top, drop = FALSE]))
    ggplot(dd, aes(exposure, signature)) +
      geom_boxplot(outlier.shape = NA, colour = "grey35", fill = NA, width = .6) +
      geom_point(position = position_jitter(height = .15, seed = 7), size = 1.6,
                 colour = TEAL, alpha = .7) +
      labs(x = "share of a cell's mutations", y = NULL,
           title = sprintf("%d cells with a signature fit", nrow(m))) + theme_lab()
  })

  output$qc_head <- renderUI({
    d <- cm()
    a <- d$ado[is.finite(d$ado)]; n <- d$n_snv[is.finite(d$n_snv)]
    bits <- sprintf("%d cells", nrow(d))
    if (length(a)) bits <- c(bits, sprintf("median allelic dropout %.1f%%", median(a)))
    if (length(n)) bits <- c(bits, sprintf("median %s SNVs called per cell",
                                           format(round(median(n)), big.mark = ",")))
    ts <- table(d$top_signature[!is.na(d$top_signature)])
    if (length(ts)) bits <- c(bits, sprintf("%s dominates in %d of %d cells",
                                            names(which.max(ts)), max(ts), sum(ts)))
    note(paste0(paste(bits, collapse = "; "), "."))
  })

  output$af_plot <- renderPlot({
    d <- AF
    d$dir <- ifelse(d$delta > 0, "rose", ifelse(d$delta < 0, "fell", "unchanged"))
    lab <- d[order(-abs(d$delta)), ][1:8, ]
    lim <- c(0, max(d$pre, d$post) * 1.08)
    ggplot(d, aes(pre, post)) +
      geom_abline(slope = 1, intercept = 0, colour = "#C7CDD2", linetype = 2) +
      geom_point(aes(colour = dir, size = cells_pre + cells_post), alpha = .9) +
      geom_text(data = lab, aes(label = gene), vjust = -1.1, size = 3.4, colour = "#17212B") +
      scale_colour_manual(values = c(rose = CARD, fell = TEAL, unchanged = GREY), name = NULL) +
      scale_size_continuous(range = c(2, 6), guide = "none") +
      coord_equal(xlim = lim, ylim = lim) +
      labs(x = "allele frequency before treatment (sample 4272)",
           y = "allele frequency after treatment (sample 4295)") + theme_lab()
  })

  output$af_tbl <- renderDT({
    d <- AF[order(-AF$delta), c("gene", "locus", "pre", "post", "delta", "cells_pre", "cells_post")]
    names(d) <- c("Gene", "Locus", "Before", "After", "Change", "Cells before", "Cells after")
    datatable(d, rownames = FALSE, options = list(pageLength = 8, dom = "ftip")) |>
      formatRound(c("Before", "After", "Change"), 3)
  })

  output$af_head <- renderUI({
    up <- AF[which.max(AF$delta), ]; dn <- AF[which.min(AF$delta), ]
    note(sprintf(paste("%d of the panel's variants have coverage in both samples. %s rises most",
                       "(%.3f to %.3f) and %s falls most (%.3f to %.3f); %d of %d rose."),
                 nrow(AF), up$gene, up$pre, up$post, dn$gene, dn$pre, dn$post,
                 sum(AF$delta > 0), nrow(AF)))
  })

  # ---- relative contribution of the mutational signatures ----
  # Exposures are per-cell fractions summing to 1, so the honest aggregate is
  # weighted by each cell's mutation count, not a mean over cells.
  sig_share <- function(patients) {
    k <- intersect(CELLMETA$cell[CELLMETA$patient %in% patients & is.finite(CELLMETA$n_snv)],
                   rownames(SIG))
    if (!length(k)) return(NULL)
    w <- CELLMETA$n_snv[match(k, CELLMETA$cell)]
    v <- colSums(SIG[k, , drop = FALSE] * w) / sum(w)
    list(share = sort(v[v > 0], decreasing = TRUE), n_cells = length(k), n_mut = sum(w))
  }

  pool <- function(v, floor = 0.02) {
    big <- v[v >= floor]
    if (sum(v < floor) > 0) big <- c(big, `other signatures` = sum(v[v < floor]))
    big
  }

  output$pie_title <- renderText(
    sprintf("Relative contribution of the mutational signatures \u2014 patient %s", input$cpat))

  output$sig_pie <- renderPlot({
    r <- sig_share(input$cpat)
    validate(need(!is.null(r), sprintf("No signature fit for patient %s.", input$cpat)))
    v <- pool(r$share, floor = 0.03)
    d <- data.frame(sig = factor(names(v), levels = names(v)), share = as.numeric(v))
    d$lab <- ifelse(d$share >= .04, sprintf("%s\n%.0f%%", d$sig, 100 * d$share), "")
    d$pos <- cumsum(d$share) - d$share / 2
    ggplot(d, aes(x = "", y = share, fill = sig)) +
      geom_col(width = 1, colour = "white", linewidth = .6) +
      coord_polar(theta = "y", direction = -1) +
      geom_text(aes(y = 1 - pos, label = lab), size = 3.4, colour = "white", fontface = "bold") +
      scale_fill_manual(values = grDevices::hcl.colors(nrow(d), "Spectral"), name = NULL) +
      theme_void(base_size = 12) + theme(legend.position = "right")
  })

  output$sig_bar <- renderPlot({
    pats <- setdiff(sort(unique(CELLMETA$patient)), "Invitro")
    rows <- lapply(pats, function(q) {
      r <- sig_share(q); if (is.null(r)) return(NULL)
      v <- pool(r$share)
      data.frame(patient = sprintf("%s\n(%s mutations)", q, format(r$n_mut, big.mark = ",")),
                 sig = names(v), share = as.numeric(v), stringsAsFactors = FALSE)
    })
    d <- do.call(rbind, rows)
    validate(need(!is.null(d) && nrow(d) > 0, "No signature fits available."))
    ord <- names(sort(tapply(d$share, d$sig, sum), decreasing = TRUE))
    ord <- c(setdiff(ord, "other signatures"), intersect("other signatures", ord))
    d$sig <- factor(d$sig, levels = rev(ord))
    ggplot(d, aes(share, patient, fill = sig)) +
      geom_col(colour = "white", linewidth = .3) +
      scale_fill_manual(values = rev(grDevices::hcl.colors(length(ord), "Spectral")), name = NULL) +
      guides(fill = guide_legend(reverse = TRUE)) +
      scale_x_continuous(labels = function(x) paste0(x * 100, "%"),
                         expand = expansion(mult = c(0, .01))) +
      labs(x = "share of that patient's mutations", y = NULL,
           title = "The four patients side by side") + theme_lab()
  })

  output$pie_head <- renderUI({
    r <- sig_share(input$cpat)
    if (is.null(r)) return(NULL)
    top <- utils::head(r$share, 3)
    note(sprintf("%s mutations pooled from %d cells. %s.", format(r$n_mut, big.mark = ","), r$n_cells,
                 paste(sprintf("%s accounts for %.0f%% of them", names(top), 100 * top), collapse = ", ")))
  })

  # ---- Phylogeny ----
  cur_tree <- reactive({
    txt <- TREES[[input$tree]]; req(!is.null(txt))
    t <- ape::read.tree(text = txt); req(!is.null(t))
    t
  })

  # Only the Figure 7 patient has timepoints, clones and copy-number calls. The
  # other trees would colour uniformly grey, so the picker offers each tree the
  # annotations that tree actually has.
  tip_keys <- reactive(norm_tip(cur_tree()$tip.label))

  avail_cols <- reactive({
    k <- tip_keys()
    opts <- c("None" = "none")
    if (any(k %in% CELLMETA$cell[!is.na(CELLMETA$top_signature)]))
      opts <- c(opts, "Signature composition (pie per cell)" = "pie",
                      "Dominant mutational signature" = "signature")
    if (any(k %in% CELLMETA$cell[is.finite(CELLMETA$ado)]))
      opts <- c(opts, "Allelic dropout rate" = "ado", "Sequencing depth" = "depth")
    if (any(k %in% TIPS$cell)) {
      if (length(unique(TIPS$timepoint[match(k, TIPS$cell)])) > 1)
        opts <- c(opts, "Timepoint" = "timepoint")
      opts <- c(opts, "Clone" = "clone", "Chromosome 4 deletion" = "chr4")
    }
    opts
  })

  output$tipcol_ui <- renderUI({
    o <- avail_cols()
    radioButtons("tipcol", "Colour tips by", choices = o,
                 selected = if (length(o) > 1) o[[2]] else o[[1]])
  })

  # one row per tip, columns pooled to the signatures that actually show at this
  # scale; a tip with no fit gets an all-zero row and is drawn as a grey dot
  tip_pie <- reactive({
    k <- tip_keys()
    m <- matrix(0, nrow = length(k), ncol = ncol(SIG), dimnames = list(k, colnames(SIG)))
    i <- match(k, rownames(SIG)); ok <- !is.na(i)
    m[ok, ] <- SIG[i[ok], , drop = FALSE]
    tot <- colSums(m)
    keep <- names(sort(tot[tot > 0], decreasing = TRUE))
    keep <- keep[seq_len(min(7L, length(keep)))]
    rest <- setdiff(colnames(m)[colSums(m) > 0], keep)
    out <- m[, keep, drop = FALSE]
    if (length(rest)) out <- cbind(out, `other` = rowSums(m[, rest, drop = FALSE]))
    rs <- rowSums(out)
    out[rs > 0, ] <- out[rs > 0, , drop = FALSE] / rs[rs > 0]   # each pie sums to 1
    list(m = out, has = rs > 0)
  })

  # Timepoint of each tip, where the tree spans both samples. Only the clone tree
  # does; the CellPhy trees are single-sample, so this is NULL for them.
  tip_timepoint <- reactive({
    v <- TIPS$timepoint[match(tip_keys(), TIPS$cell)]
    if (length(unique(stats::na.omit(v))) > 1) as.character(v) else NULL
  })

  TP_EDGE <- c(Pretreatment = "#4A555F", `Post-treatment` = CARD)

  BINS <- function(x, n = 5) {
    ok <- is.finite(x); if (!any(ok)) return(rep(NA_character_, length(x)))
    b <- unique(stats::quantile(x[ok], seq(0, 1, length.out = n + 1), na.rm = TRUE))
    if (length(b) < 3) return(ifelse(ok, sprintf("%.2f", x), NA))
    as.character(cut(x, breaks = b, include.lowest = TRUE, dig.lab = 3))
  }

  tip_ann <- reactive({
    k <- tip_keys(); col <- input$tipcol
    req(!is.null(col))
    j <- match(k, CELLMETA$cell); i <- match(k, TIPS$cell)
    v <- switch(col,
      signature = CELLMETA$top_signature[j],
      ado       = BINS(CELLMETA$ado[j]),
      depth     = BINS(CELLMETA$depth[j]),
      timepoint = TIPS$timepoint[i],
      clone     = TIPS$clone[i],
      chr4      = ifelse(is.na(TIPS$chr4[i]), NA,
                         ifelse(TIPS$chr4[i] == 1, "deleted", "retained")),
      rep("cell", length(k)))
    ifelse(is.na(v), "not annotated", as.character(v))
  })

  output$tree_title <- renderText({
    t <- cur_tree()
    sprintf("%s: %d cells", names(TREE_CHOICES)[match(input$tree, TREE_CHOICES)], ape::Ntip(t))
  })

  # The legend is HTML, not part of the graphic. ape::plot.phylo manages its own
  # margins, so anything drawn with legend() lands on the caption below the image.
  tree_pal <- reactive({
    if (identical(input$tipcol, "pie")) {
      pm <- tip_pie()$m
      setNames(grDevices::hcl.colors(ncol(pm), "Spectral"), colnames(pm))
    } else {
      lv <- sort(unique(tip_ann()))
      if (identical(input$tipcol, "timepoint"))
        c(Pretreatment = "#4A555F", `Post-treatment` = CARD, `not annotated` = GREY)[lv]
      else if (input$tipcol %in% c("ado", "depth")) {
        o <- lv[order(suppressWarnings(as.numeric(sub("^[\\[(]", "", sub(",.*", "", lv)))))]
        setNames(grDevices::hcl.colors(length(o), "Viridis"), o)[lv]
      } else setNames(grDevices::hcl.colors(length(lv), "Spectral"), lv)
    }
  })

  output$tree_legend <- renderUI({
    dot <- function(col) span(style = sprintf(
      "display:inline-block;width:.75rem;height:.75rem;border-radius:50%%;background:%s;margin-right:.3rem;vertical-align:-1px", col))
    bar <- function(col) span(style = sprintf(
      "display:inline-block;width:1.1rem;height:.2rem;background:%s;margin-right:.3rem;vertical-align:.18rem", col))
    rows <- list()
    if (!identical(input$tipcol, "none")) {
      pal <- tree_pal(); pal[is.na(pal)] <- GREY
      rows[[length(rows) + 1]] <- div(
        style = "display:flex;flex-wrap:wrap;gap:.35rem 1rem;margin:.4rem 0 0;font-size:.85rem",
        span(style = "color:#64707C", "tips:"),
        lapply(names(pal), function(n) span(dot(pal[[n]]), n)))
    }
    if (!is.null(tip_timepoint()))
      rows[[length(rows) + 1]] <- div(
        style = "display:flex;flex-wrap:wrap;gap:.35rem 1rem;margin:.3rem 0 0;font-size:.85rem",
        span(style = "color:#64707C", "terminal branch:"),
        lapply(names(TP_EDGE), function(n) span(bar(TP_EDGE[[n]]), n)),
        span(bar("#C7CDD2"), "sample unknown"))
    do.call(tagList, rows)
  })

  output$tree_plot <- renderPlot(height = function() {
    n <- tryCatch(ape::Ntip(cur_tree()), error = function(e) 40)
    per <- if (identical(input$tipcol, "pie")) 30 else 9   # pies need room to read
    max(560, min(3600, round(per * n)))
  }, {
    t <- cur_tree(); ann <- tip_ann()
    lv <- sort(unique(ann))
    pal <- if (identical(input$tipcol, "timepoint"))
             c(Pretreatment = "#4A555F", `Post-treatment` = CARD, `not annotated` = GREY)[lv]
           else if (input$tipcol %in% c("ado", "depth")) {
             o <- lv[order(suppressWarnings(as.numeric(sub("^[\\[(]", "", sub(",.*", "", lv)))))]
             setNames(grDevices::hcl.colors(length(o), "Viridis"), o)[lv]
           } else setNames(grDevices::hcl.colors(length(lv), "Spectral"), lv)
    pal[is.na(pal)] <- GREY
    cols <- pal[ann]
    # The terminal branch of each cell is coloured by the sample that cell came
    # from, so before and after treatment are readable straight off the topology.
    tp <- tip_timepoint()
    ecol <- rep("#5A6570", nrow(t$edge)); ewid <- rep(.9, nrow(t$edge))
    if (!is.null(tp)) {
      term <- t$edge[, 2] <= ape::Ntip(t)
      v <- tp[t$edge[term, 2]]
      ecol[term] <- ifelse(is.na(v), "#C7CDD2", unname(TP_EDGE[v]))
      ewid[term] <- 2.1
    }
    par(mar = c(1, 1, 1, 1), xpd = TRUE)
    plot(t, type = input$ttype, show.tip.label = isTRUE(input$tiplab),
         tip.color = cols, cex = .6, no.margin = FALSE,
         edge.color = ecol, edge.width = ewid,
         use.edge.length = !is.null(t$edge.length))
    if (identical(input$tipcol, "pie")) {
      tp <- tip_pie(); pm <- tp$m
      pcol <- unname(tree_pal())
      # scale the pies to the tip count so a 113-tip tree does not turn to soup
      # A pie's radius scales with the plot WIDTH, not the tip count, so the size is
      # a constant wherever each tip gets its intended 30 px of height. Scaling it
      # by 1/n, as a first attempt did, shrank an 85-tip tree's pies for no reason.
      # Where the height cap bites, fall back to the room a tip actually has.
      per_tip <- (par("din")[2] * 96) / max(ape::Ntip(t), 1)
      pc <- 0.70 * min(1, per_tip / 30)
      if (any(!tp$has)) ape::tiplabels(pch = 19, col = GREY, cex = .6,
                                       tip = which(!tp$has))
      if (any(tp$has)) ape::tiplabels(pie = pm[tp$has, , drop = FALSE],
                                      tip = which(tp$has), piecol = pcol, cex = pc)
    } else {
      ape::tiplabels(pch = 19, col = cols, cex = 1.1)
    }
    # bootstrap support, only where the file recorded it and it clears the threshold
    if (!is.null(t$node.label)) {
      b <- suppressWarnings(as.numeric(t$node.label))
      keep <- is.finite(b) & b >= input$bootmin
      if (any(keep)) ape::nodelabels(text = as.character(b[keep]), node = which(keep) + ape::Ntip(t),
                                     frame = "none", col = "#17212B", cex = .6, adj = c(1.1, -.3))
    }
  })

  output$tree_head <- renderUI({
    t <- cur_tree(); ann <- tip_ann()
    bits <- sprintf("%d cells", ape::Ntip(t))
    if (!is.null(t$edge.length))
      bits <- c(bits, sprintf("total branch length %.3f substitutions per site", sum(t$edge.length)))
    else bits <- c(bits, "no branch lengths in this tree, so it is drawn as a cladogram")
    if (!is.null(t$node.label)) {
      b <- suppressWarnings(as.numeric(t$node.label)); b <- b[is.finite(b)]
      if (length(b)) bits <- c(bits, sprintf("%d of %d internal nodes reach %d%% bootstrap support",
                                             sum(b >= input$bootmin), length(b), input$bootmin))
    }
    na <- if (identical(input$tipcol, "pie")) sum(!tip_pie()$has) else sum(ann == "not annotated")
    if (na) bits <- c(bits, sprintf("%d %s no annotation for this colouring",
                                    na, if (na == 1) "tip has" else "tips have"))
    note(paste0(paste(bits, collapse = "; "), "."))
  })

  # ---- Genes ----
  gene_am <- reactive({
    a <- stats::aggregate(am_pathogenicity ~ gene, am, mean)
    names(a)[2] <- "mean_am"
    m <- merge(rec, a, by = "gene", all.x = TRUE)
    m$label <- GRP_LAB[m$group]
    m
  })

  output$gene_scatter <- renderPlot({
    d <- gene_am(); d <- d[is.finite(d$mean_am) & is.finite(d$fold_size_corr), ]
    sel <- d[d$gene == input$gene, ]
    ggplot(d, aes(fold_size_corr, mean_am, colour = label)) +
      geom_hline(yintercept = c(.34, .564), linetype = "22", colour = "grey55", linewidth = .4) +
      geom_point(size = 3, alpha = .9) +
      scale_x_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300, 700)) +
      { if (nrow(sel)) geom_point(data = sel, size = 6, shape = 21, fill = NA,
                                  colour = "black", stroke = 1.1) } +
      { if (nrow(sel)) geom_text(data = sel, aes(label = gene), vjust = -1.4,
                                 fontface = "bold", show.legend = FALSE) } +
      scale_colour_manual(values = setNames(GRP[names(GRP_LAB)], GRP_LAB), name = NULL) +
      labs(x = "recurrence, fold above the coding-size expectation",
           y = "mean AlphaMissense pathogenicity") + theme_lab()
  })

  output$gene_count <- renderUI({
    d <- gene_am()
    n <- sum(is.finite(d$mean_am) & is.finite(d$fold_size_corr))
    miss <- d$gene[!(is.finite(d$mean_am) & is.finite(d$fold_size_corr))]
    note(sprintf("%d of %d genes can be placed. %s carry no scored missense variant or no coding length, so they have no position here.",
                 n, nrow(d), paste(sort(miss), collapse = ", ")))
  })

  output$am_title <- renderText(sprintf("%s: every scored missense variant", input$gene))
  output$am_plot <- renderPlot({
    d <- am[am$gene == input$gene, ]
    validate(need(nrow(d) > 0, sprintf("No scored missense variants observed in %s.", input$gene)))
    d$tp <- factor(ifelse(d$timepoint == "relapse", "Relapse", "Diagnosis"),
                   levels = c("Diagnosis", "Relapse"))
    ggplot(d, aes(tp, am_pathogenicity)) +
      geom_hline(yintercept = c(.34, .564), linetype = "22", colour = "grey55", linewidth = .4) +
      geom_boxplot(width = .5, outlier.shape = NA, colour = "grey30", fill = NA) +
      geom_point(position = position_jitter(width = .12, seed = 4), size = 2.4,
                 shape = 21, fill = "white", colour = "grey30") +
      ylim(0, 1.02) + labs(x = NULL, y = "AlphaMissense pathogenicity") + theme_lab()
  })

  output$gene_tbl <- renderDT({
    d <- gene_am()[, c("gene", "label", "nonsyn", "fold_size_corr", "relapse_freq_pct", "mean_am")]
    names(d) <- c("Gene", "Class", "Nonsyn. events", "Fold recurrence", "Relapse %", "Mean AlphaMissense")
    datatable(d[order(-d$`Fold recurrence`), ], rownames = FALSE, selection = "single",
              options = list(pageLength = 8, dom = "ftip")) |>
      formatRound(c("Fold recurrence", "Relapse %", "Mean AlphaMissense"), 2)
  })
  observeEvent(input$gene_tbl_rows_selected, {
    i <- input$gene_tbl_rows_selected
    if (length(i)) {
      d <- gene_am(); d <- d[order(-d$fold_size_corr), ]
      updateSelectizeInput(session, "gene", selected = d$gene[i])
    }
  })
}

shinyApp(ui, server)
