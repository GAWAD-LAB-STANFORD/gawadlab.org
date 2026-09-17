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
BULK <- B$bulk; WGSMUT <- B$wgs

# The four patients whose Figure 5 tree mutations were re-measured in bulk at
# both timepoints, and whose trees carry exome-mapped mutations on their branches.
WGS_PT <- names(WGSMUT)

# Put a branch annotation back on a branch. The annotation was keyed on the tip
# set beneath the branch in the TreeMut object, which carries an extra all-
# reference tip ("zeros") and is rooted differently from the CellPhy tree drawn
# here. Both describe the same unrooted topology, so match on the SPLIT instead:
# drop "zeros", and name the split by its smaller side so that naming a clade or
# its complement gives the same key. Verified to place all 21 annotations across
# the four trees, with no split occurring twice.
node_sets <- function(t) {
  nt <- ape::Ntip(t); kids <- split(t$edge[, 2], t$edge[, 1])
  out <- vector("list", nt + ape::Nnode(t))
  for (i in seq_len(nt)) out[[i]] <- t$tip.label[i]
  for (n in rev(sort(unique(t$edge[, 1]))))
    out[[n]] <- sort(unlist(out[as.integer(kids[[as.character(n)]])]))
  out
}
split_key <- function(s, all) {
  s <- setdiff(s, "zeros"); o <- setdiff(all, s)
  paste(if (length(s) <= length(o)) s else o, collapse = "|")
}
CELLMETA <- B$cellmeta; SIG <- B$sig; POS <- B$pos

# Newick is parsed here rather than shipped as a serialised tree, so the bundle
# stays in kilobytes. Labels come through as 4295_A10 or 4295.F1 depending on the
# file; the annotation table is keyed on the normalised 4295-A10 form.
TREE_CHOICES <- local({
  n <- vapply(TREES, function(x) length(gregexpr("[,(]", x)[[1]]) , integer(1))
  lab <- c("4295" = "Patient 4295", "445" = "Patient 445", "417" = "Patient 417",
           "4084" = "Patient 4084", "Invitro" = "In vitro benchmark",
           "clone" = "Patient 4295, before and after induction")
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

# Lay out labels in a right-hand column without overlaps. ggrepel is not in the
# webR repo, so: push each label below the one above by a fixed gap, then move
# the whole stack to fit between lo and hi - translating it, or compressing it
# when there are more labels than the panel has room for. Never clamp: clamping
# is what stacks several labels on the same limit. Returns values in the order
# given, so the caller need not pre-sort.
stack_labels <- function(v, lo, hi, frac = .04) {
  n <- length(v); if (!n) return(numeric(0))
  span <- hi - lo
  if (!is.finite(span) || span <= 0) return(v)
  gap <- if (n > 1) min(span * frac, span / (n - 1)) else 0
  o <- order(-v); y <- v[o]
  for (i in seq_len(n)[-1]) if (y[i - 1] - y[i] < gap) y[i] <- y[i - 1] - gap
  r <- range(y)
  if (diff(r) > span)      y <- lo + (y - r[1]) / diff(r) * span
  else if (r[1] < lo)      y <- y + (lo - r[1])
  else if (r[2] > hi)      y <- y - (r[2] - hi)
  out <- numeric(n); out[o] <- y; out
}
ext  <- function(l, u) a(l, href = u, target = "_blank", rel = "noopener")

# A browser is only useful to a researcher if the numbers can leave it, so every
# table ships the exact rows and columns on screen as a CSV. DT's own Buttons
# extension is not bundled by shinylive, so this goes through downloadHandler,
# which runs in the wasm VM and hands the browser a real file.
dl_link <- function(id) div(
  style = "margin-top:.45rem",
  downloadLink(paste0("dl_", id), "Download this table (CSV)",
               style = "font-size:.85rem;color:#64707C"))


# ---------------------------------------------------------------------------
# The same patients appear under three different identifiers in this study: the
# short bulk name on the drug-response experiment (1678), the St Jude accession
# on the RAS and targeted work (SJETV022), and the numbered patient in the
# manuscript (patient 14). Table S1 of the paper is the cross-reference; it is
# reproduced here so a reader never has to guess whether two tabs are showing
# the same person. Drug-response 3072 and RAS SJETV026 are one patient.
# ---------------------------------------------------------------------------
PTKEY <- data.frame(
  id = c("1678","2364","2488","2788","3072",
         "SJETV022","SJETV024","SJETV025","SJETV026","SJETV075","SJETV077","SJETV078",
         "SJETV083","SJETV092","4295","417","445","4084"),
  sj = c("SJETV022","SJETV024","SJETV078","SJETV025","SJETV026",
         "SJETV022","SJETV024","SJETV025","SJETV026","SJETV075","SJETV077","SJETV078",
         "SJETV083","SJETV092","-","-","-","-"),
  bulk = c("1678","2364","2488","2788","3072",
           "1678","2364","2788","3072","2185","2295","2488","1178","0060","-","-","-","-"),
  patient = c(14,15,5,8,2, 14,15,8,2,3,4,5, 6,7, 16,17,18,19),
  subtype = c(rep("ETV6-RUNX1", 14),
              "Ph-like (CRLF2, JAK-mutant)","Ph-like (IGH-CRLF2)",
              "Hypodiploid (45, -X, -7)","No recurrent lesion (normal karyotype)"),
  stringsAsFactors = FALSE)
pt_label <- function(x) {
  i <- match(as.character(x), PTKEY$id)
  ifelse(is.na(i), as.character(x),
         ifelse(PTKEY$sj[i] == "-" | PTKEY$sj[i] == as.character(x),
                sprintf("%s  (patient %d)", x, PTKEY$patient[i]),
                sprintf("%s  = %s, patient %d", x, PTKEY$sj[i], PTKEY$patient[i])))
}
pt_choices <- function(v) { v <- sort(unique(as.character(v))); setNames(v, pt_label(v)) }

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
      selectInput("dpat", "Patient", pt_choices(drug$patient)),
      selectizeInput("dmut", "Highlight a mutation", choices = colnames(DMAT),
                     options = list(maxOptions = 200)),
      hr(),
      radioButtons("sel_drug", "Rank mutations selected by",
                   c("Daunorubicin (DNR-Hi)" = "DNR-Hi",
                     "Prednisolone (Pred-Hi)" = "Pred-Hi")),
      radioButtons("sel_scope", "Across",
                   c("All five patients" = "all", "This patient only" = "one")),
      checkboxInput("sel_consist", "Only where every treated replicate exceeds every DMSO replicate", FALSE)),
    conditionalPanel("input.nav == 'Single cells'",
      selectInput("cpat", "Patient", pt_choices(CELLMETA$patient)),
      radioButtons("cfill", "Colour the paired-sample plot by",
                   c("Clone" = "clone", "Timepoint" = "timepoint",
                     "Chromosome 4 deletion" = "Chr4_Deletion",
                     "Chromosome 6 deletion" = "Chr6_1_Deletion"))),
    conditionalPanel("input.nav == 'Induction therapy'",
      radioButtons("ind_mark", "Tip marks",
                   c("Mutational signature (cells after induction only)" = "sig",
                     "Sample the cell came from" = "tp"),
                   selected = "sig"),
      radioButtons("ind_type", "Layout",
                   c("Cladogram" = "phylogram", "Fan" = "fan", "Unrooted" = "unrooted")),
      checkboxInput("ind_node", "Clade pies: pre/post composition at each node", TRUE),
      sliderInput("ind_min", "Smallest clade to draw a pie for", 2, 40, 6, step = 1),
      checkboxInput("ind_lab", "Show cell labels", FALSE)),
    conditionalPanel("input.nav == 'Phylogeny'",
      selectInput("tree", "Tree", TREE_CHOICES),
      radioButtons("ttype", "Layout",
                   c("Phylogram" = "phylogram", "Fan" = "fan", "Unrooted" = "unrooted")),
      uiOutput("tipcol_ui"),
      sliderInput("bootmin", "Label bootstrap support of at least", 0, 100, 50, step = 5),
      checkboxInput("treemut", "Mark branches carrying mutations", TRUE),
      checkboxInput("tiplab", "Show cell labels", FALSE)),
    conditionalPanel("input.nav == 'Genes'",
      selectizeInput("gene", "Highlight a gene", choices = sort(unique(rec$gene)),
                     selected = if ("TBL1XR1" %in% rec$gene) "TBL1XR1" else sort(unique(rec$gene))[1],
                     options = list(maxOptions = 200))),
    hr(),
    note(strong("Pang, Prieto ", em("et al."), "."),
         " Single-cell sequencing reveals extensive genetic diversity underlying ",
         "pediatric ALL treatment complexity. Every panel is drawn from the ",
         "manuscript's figure-level source data; nothing here is simulated. ",
         "Manuscript in review, so there is no DOI to cite yet.",
         div(style = "margin-top:.5rem;display:flex;flex-direction:column;gap:.15rem",
             ext("Source data and count matrices", "https://gawadlab.org/apps.html"),
             ext("Code on GitHub", "https://github.com/GAWAD-LAB-STANFORD"),
             ext("Gawad Lab", "https://gawadlab.org")),
         div(style = "margin-top:.5rem",
             strong("Reading the sample names. "),
             "Each patient was sampled twice, and is named for the later sample. ",
             "For patient 4295 the pair is 4272 before induction and 4295 after it; ",
             "the others are 368/417, 380/445 and 4072/4084, earlier number first. ",
             "Before and after induction mean the same thing on every tab: before is ",
             "diagnosis, after is four weeks of induction therapy."),
         div(style = "margin-top:.5rem",
             "Every table on every tab has a CSV download beneath it."))
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
           note("SNVs per megabase in a published survey of 961 paediatric tumours. ",
                "ALL sits at the bottom of the range, which is the observation this ",
                "study set out to reconcile with the complexity of its treatment. ",
                "The axis is logarithmic, so the 12 tumours reported at zero SNVs per ",
                "megabase cannot be drawn and are absent from the plot."))),
    card(card_header("Per-patient detail"), tagList(DTOutput("burden_tbl"), dl_link("burden_tbl")))),

  nav_panel("RAS",
    card(card_header("Activating RAS mutations found by error-corrected sequencing"),
         plotOutput("ras_plot", height = 420),
         note("Each point is one activating mutation. Bulk sequencing reported a single ",
              "RAS mutation in each of these patients; error-corrected sequencing finds ",
              "several more at low allele frequency, which is why they were missed.")),
    layout_columns(col_widths = c(7, 5),
      card(card_header("Allele frequency by codon"), plotOutput("ras_codon", height = 340)),
      card(card_header("Every RAS mutation"), tagList(DTOutput("ras_tbl"), dl_link("ras_tbl"))))),

  nav_panel("Drug response",
    card(card_header(textOutput("sel_plot_title")),
         plotOutput("sel_plot", height = 520),
         uiOutput("sel_note")),

    card(card_header(textOutput("drug_title")), plotOutput("drug_heat", height = 460),
         note("Mutant allele frequency in percent, one row per sequenced sample and one ",
              "column per recurrent mutation. Pred-Hi is prednisolone and DNR-Hi is ",
              "daunorubicin, each against its own DMSO control and the diagnostic bulk ",
              "sample. A column that rises under one drug and not the other marks a ",
              "population with differential sensitivity. ", tags$b("Grey tiles are not missing data"),
              " - every tile was measured. Grey means the variant was not called in that sample; "
              , "called frequencies in this experiment start at 15%, so grey reads as ",
              tags$em("below the calling threshold"), " rather than as proven absent.")),
    layout_columns(col_widths = c(6, 6),
      card(card_header(textOutput("dmut_title")),
           note(tags$b("All five patients, side by side,"), " so one mutation can be compared ",
                "across the cohort; the selected patient's panel is boxed in amber."),
           plotOutput("dmut_plot", height = 340)),
      card(card_header("Where this mutation sits in the composite ranking"),
           uiOutput("dmut_rank_note"))),

    card(card_header("Patient key - the same patients appear under three identifiers"),
         DTOutput("ptkey_tbl"),
         note("From Table S1 of the manuscript. The drug-response experiment names ",
              "patients by the short bulk sample name, the RAS and targeted work uses the ",
              "St Jude accession, and the manuscript numbers them. Drug-response ",
              tags$b("3072"), " and RAS ", tags$b("SJETV026"), " are the same patient, as are ",
              tags$b("2488"), " and ", tags$b("SJETV078"), ", and ", tags$b("2788"),
              " and ", tags$b("SJETV025"), ".")),
    card(card_header(textOutput("sel_title")),
         tagList(DTOutput("sel_tbl"), dl_link("sel_tbl")))),

  nav_panel("SJETV077 panel",
    card(card_header("SJETV077 = sample 2295, patient 4 - across nine ex vivo conditions"),
         plotOutput("sj_plot", height = 620),
         note("A separate single-patient experiment covering six agents plus controls. ",
              "It has its own tab because it is not one of the five patients on the Drug ",
              "response tab, so no patient selector applies to it. Grey tiles are ",
              "combinations with no reported measurement, not zeros - roughly half of this ",
              "grid was not reported."))),

  nav_panel("Single cells",
    card(card_header(textOutput("qc_title")),
      layout_columns(col_widths = c(6, 6),
        plotOutput("qc_plot", height = 380), plotOutput("sig_plot", height = 380)),
      uiOutput("qc_head"),
      note("Every patient's cells are here. Allelic dropout and depth are the two ",
           "measurements that decide whether a single-cell genome can be called at all, ",
           "and the signature panel is the COSMIC exposure fitted to each cell's own ",
           "mutations. The panels below need the paired before- and after-induction ",
           "samples, which exist for one patient only.")),

    layout_columns(col_widths = c(7, 5),
      card(card_header(textOutput("cell_plot_title")),
           uiOutput("cell_plot_slot"),
           note("Two samples from patient 4295: 4272 drawn before induction (30 cells) and ",
                "4295 drawn after it (85 cells). Axes are the measured surface-marker intensities used to ",
                "separate leukemic from normal and premalignant cells.")),
      card(card_header(textOutput("clone_plot_title")), uiOutput("clone_plot_slot"),
           note("Each bar is the percentage of that sample's cells, not a raw count, ",
                "because 30 cells were sequenced before induction against 85 after. Clones ",
                strong(paste(EMERGENT, collapse = " and ")),
                " are absent from the before-induction sample entirely and appear only after ",
                "it, but they are 3 and 1 cells, so read them as the observation ",
                "they are rather than as a reliable frequency."))),

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
         uiOutput("tree_legend"), uiOutput("tree_muts"), uiOutput("tree_head"),
         note("The four patient trees hold cells from after induction only: every cell in them was ",
              "taken after four weeks of induction therapy, which is why their branches are ",
              "uniform. Only the last tree spans two samples from patient 4295, 4272 drawn ",
              "before induction (29 cells) and 4295 after (84 cells), and there the branch into ",
              "each cell carries the sample it came from. A tip pie is that one cell's own ",
              "mutations split by signature. ",
              "The patient trees are maximum-likelihood phylogenies built by CellPhy from ",
              "somatic single-nucleotide variants, with support from 100 bootstrap replicates; ",
              "branch lengths are substitutions per site. The last tree is the topology ",
              "behind Figure 7A, which carries no branch lengths, so it is drawn as ",
              "a cladogram and the fan and unrooted layouts show topology only."))),

  nav_panel("Induction therapy",
    card(card_header("Two strategies, two kinds of before and after"),
      note("Induction therapy appears twice in this paper, measured two different ways. ",
           "Figure 5 sequenced whole genomes of single cells with PTA and built four ",
           "maximum-likelihood trees with CellPhy \u2014 but every one of those 165 cells was ",
           "drawn ", tags$b("after"), " four weeks of induction, from patients who still had ",
           "detectable residual disease. Those trees have no before-and-after inside them. ",
           "Their before-and-after is bulk: each tree's mutations were re-measured in the ",
           "diagnostic sample and again in the remission sample, which is the next two panels ",
           "here, and it covers all four patients. Figure 7 did something different \u2014 ",
           "single-cell ", tags$b("exome"), " sequencing with a ConDoR tree \u2014 and only there ",
           "do single cells from before and after sit in the same phylogeny. That is one ",
           "patient, 4295, and it is the rest of this tab. The trees themselves are on the ",
           "Phylogeny tab.")),

    card(card_header(textOutput("wgs_card_title")),
      radioButtons("wgs_pt", NULL, inline = TRUE,
                   choices = setNames(WGS_PT, pt_label(WGS_PT)),
                   selected = if ("4295" %in% WGS_PT) "4295" else WGS_PT[1]),
      plotOutput("wgs_slope", height = 480),
      uiOutput("wgs_slope_head"),
      note("Every mutation that could be placed on a branch of this patient's tree on the ",
           "previous tab, measured in the diagnostic bulk ",
           "sample and again in the remission bulk sample. The y axis is square-root scaled so ",
           "that the small surviving frequencies stay visible next to the clonal ones. A line ",
           "falling to the floor is a mutation that became undetectable in the bulk \u2014 which ",
           "is not the same as gone, because a bulk remission sample is mostly normal marrow ",
           "and the residual cells are rare. That gap is exactly why the single-cell work in ",
           "Figure 7 was needed.")),

    card(card_header(textOutput("wgs_tbl_title")), tagList(DTOutput("wgs_tbl"), dl_link("wgs_tbl"))),

    card(class = "border-0 bg-transparent",
         note(tags$b("Everything below is the single-cell exome experiment."), " It was run ",
              "on patient 4295 only - the one patient sequenced cell by cell both before and ",
              "after induction - so selecting any other patient above turns these panels into ",
              "a note saying so rather than leaving 4295 on screen.")),

    card(fill = FALSE, card_header(textOutput("ind_title")),
         uiOutput("ind_tree_slot"),
         uiOutput("ind_legend_slot"), uiOutput("ind_head_slot"),
         uiOutput("ind_desc_slot")),

    card(card_header(textOutput("clade_hdr")),
         uiOutput("ind_pies_slot"), uiOutput("ind_pies_head_slot"),
         note("One pie per clade, numbered as in Figure 7A, which uses the same internal-node ",
              "indices. A pie that is entirely pink is a clade whose cells were all found after ",
              "induction; entirely blue means the clade did not survive it. Use the slider to ",
              "set how small a clade still earns a pie.")),

    card(card_header(textOutput("pospct_hdr")),
      uiOutput("pos_plot_slot"),
      uiOutput("pos_head_slot"),
      note("Each line is one variant, and the axis is the percentage of that ",
           "timepoint's cells that carry it \u2014 a carrier frequency, not an allele ",
           "frequency. Cells are all 115 with a genotype call, 30 before induction and ",
           "85 after; calls come from the same ConDoR genotype matrix Figure 7 is drawn ",
           "from, where every one of these variants is called present or absent. Two of ",
           "those cells are not tips of the Figure 7 tree, which is why the tree on the ",
           "previous tab holds 113. Variants drawn in red were carried by no cell before ",
           "induction and appear only after it.")),

    card(card_header(textOutput("postbl_hdr")), uiOutput("pos_tbl_slot")),

    card(card_header(textOutput("afplot_hdr")),
      uiOutput("af_plot_slot"),
      uiOutput("af_head_slot"),
      note("Each point is one somatic variant, pooling the alt and total reads of every cell ",
           "in that sample, so this is a pseudobulk allele frequency rather than a per-cell ",
           "call. Points above the diagonal rose under treatment. Read counts come from a ",
           "targeted panel of 31 variants, so this is not a genome-wide survey, and the ",
           "before-induction sample carries 30 cells against 85 after, which makes the ",
           "before-induction estimate the noisier of the two.")),

    card(card_header(textOutput("aftbl_hdr")), uiOutput("af_tbl_slot")),

    card(card_header(textOutput("geno_title")), uiOutput("geno_plot_slot"),
         note("Presence or absence of each somatic variant in each cell, cells ordered by ",
              "timepoint then clone. These are the called genotypes from the ConDoR matrix, ",
              "not dropout-corrected, so a blank cell means the variant was not called in ",
              "that cell rather than that it is certainly absent.")),

    card(card_header(textOutput("clone2_title")), uiOutput("clone_plot2_slot"),
         note("Each bar is the percentage of that sample's cells, not a raw count, because 30 ",
              "cells were sequenced before induction against 85 after."))
  ),

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
      card(card_header("Gene table"), tagList(DTOutput("gene_tbl"), dl_link("gene_tbl")))))
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
    # Ratio of the two pooled medians, which is NOT a paired statistic; the
    # per-patient ratio is computed separately rather than implied by wording.
    pp <- vapply(split(burden, burden$patient), function(d) {
      bb <- d$corrected_som_per_mb[d$assay == "Bulk"]
      ss <- d$corrected_som_per_mb[d$assay == "Single"]
      if (!length(bb) || !length(ss)) NA_real_ else median(ss) / median(bb)
    }, numeric(1))
    pp <- pp[is.finite(pp)]
    note(sprintf(paste("Pooling all cells and all bulk samples, the median single cell sits at",
                       "%.2f corrected mutations per Mb against %.2f for the bulks, a %.1f-fold",
                       "difference across %d cells and %d bulk samples. Taken patient by patient",
                       "instead, which is the paired comparison, the ratio runs %.1f to %.1f with",
                       "a median of %.1f across %d patients."),
                 median(s), median(b), median(s) / median(b), length(s), length(b),
                 min(pp), max(pp), median(pp), length(pp)))
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

  burden_tbl_df <- reactive({
    d <- burden[, c("patient", "assay", "total_somatic", "unique", "shared",
                    "sensitivity", "corrected_som_per_mb", "corrected_total")]
    names(d) <- c("Patient", "Assay", "Total somatic", "Unique", "Shared",
                  "Sensitivity", "Corrected som./Mb", "Corrected total")
    d
  })
  # CSV of exactly what each table shows, named so a folder of them stays legible.
  # The *_df reactives are defined further down, so the handlers resolve them by
  # name against the server environment at download time rather than now.
  SRV <- environment()
  local({
    specs <- list(
      burden_tbl = "pALL_mutation_burden",
      ras_tbl    = "pALL_RAS_mutations",
      wgs_tbl    = "pALL_branch_mutations_bulk_before_after",
      pos_tbl    = "pALL_carrier_frequency_before_after",
      af_tbl     = "pALL_allele_frequency_before_after",
      gene_tbl   = "pALL_gene_recurrence_alphamissense",
      sel_tbl    = "pALL_mutations_ranked_by_drug_selection")
    for (id in names(specs)) local({
      i <- id; stem <- specs[[id]]
      output[[paste0("dl_", i)]] <- downloadHandler(
        filename = function() sprintf("%s_%s.csv", stem, format(Sys.Date(), "%Y%m%d")),
        content  = function(file)
          utils::write.csv(get(paste0(i, "_df"), envir = SRV)(), file, row.names = FALSE))
    })
  })

  output$burden_tbl <- renderDT({
    datatable(burden_tbl_df(), rownames = FALSE, options = list(pageLength = 10, dom = "ftip")) |>
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

  ras_tbl_df <- reactive({
    d <- ras[order(-ras$AF), c("Patient", "Ras", "Location", "AA_Change", "AF")]
    names(d) <- c("Patient", "Gene", "Codon", "Change", "Allele frequency")
    d
  })
  output$ras_tbl <- renderDT({
    datatable(ras_tbl_df(), rownames = FALSE, options = list(pageLength = 8, dom = "ftip")) |>
      formatPercentage("Allele frequency", 2)
  })

  # ---- Drug response ----
  dsub <- reactive({ i <- drug$patient == input$dpat; list(meta = drug[i, ], m = DMAT[i, , drop = FALSE]) })

  output$drug_title <- renderText(sprintf("%s: ex vivo response of every recurrent mutation", pt_label(input$dpat)))
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
    # Nothing here is missing: every tile has a number. A zero means the variant was
    # not called in that sample, and the called values start at 15% - there is nothing
    # between. Showing zero as the pale end of a gradient implies a measured low
    # frequency, so it gets its own flat colour instead.
    d$af <- ifelse(d$af == 0, NA_real_, d$af)
    ggplot(d, aes(mutation, sample, fill = af)) +
      geom_tile(colour = "white", linewidth = .25) +
      scale_fill_viridis_c(option = "rocket", direction = -1, name = "AF %",
                           na.value = "#E8ECEF") +
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
      # the facet for the selected patient is boxed, so this panel still responds
      # to the Patient selector rather than sitting unchanged
      facet_wrap(~ patient, nrow = 1) +
      theme(strip.background = element_rect(
        fill = ifelse(levels(factor(d$patient)) == input$dpat, "#F6A30C", "grey92"),
        colour = NA)) +
      scale_colour_manual(values = COND_COL, guide = "none") +
      labs(x = NULL, y = "mutant allele frequency (%)") + theme_lab(45)
  })

  # Which mutations did each drug select for? The heatmap and the single-mutation
  # viewer can only answer this one mutation at a time, which is no way to search
  # 79 of them. This ranks every mutation by its rise over that patient's own DMSO
  # control, and says plainly how many such rises chance alone would produce.
  sel_rank <- reactive({
    drg <- input$sel_drug
    pats <- if (identical(input$sel_scope, "one")) input$dpat else unique(drug$patient)
    out <- list(); n_testable <- 0
    for (p in pats) {
      i <- drug$patient == p
      meta <- drug[i, ]; mm <- DMAT[i, , drop = FALSE]
      keep <- colSums(mm > 0, na.rm = TRUE) > 0
      if (!any(keep)) next
      mm <- mm[, keep, drop = FALSE]
      ctrl <- meta$condition == "DMSO"; tr <- meta$condition == drg
      if (!any(ctrl) || !any(tr)) next
      for (gmut in colnames(mm)) {
        cv <- mm[ctrl, gmut]; tv <- mm[tr, gmut]
        if (any(tv > 0, na.rm = TRUE) && length(unique(c(cv, tv))) > 1) n_testable <- n_testable + 1
        out[[length(out) + 1]] <- data.frame(
          Patient = p, Mutation = gmut,
          `DMSO %` = round(mean(cv, na.rm = TRUE), 1),
          `Treated %` = round(mean(tv, na.rm = TRUE), 1),
          `Rise` = round(mean(tv, na.rm = TRUE) - mean(cv, na.rm = TRUE), 1),
          `Every replicate above` = all(tv > max(cv, na.rm = TRUE), na.rm = TRUE),
          Flag = paste(c(if (grepl(";", gmut)) "multi-mapped",
                         if (grepl("p\\.([A-Z])[0-9]+\\1$", gmut)) "synonymous"),
                       collapse = ", "),
          check.names = FALSE, stringsAsFactors = FALSE)
      }
    }
    d <- if (length(out)) do.call(rbind, out) else
      data.frame(Patient = character(0), Mutation = character(0))
    if (nrow(d) && input$sel_consist) d <- d[d$`Every replicate above`, ]
    if (nrow(d)) d <- d[order(-d$Rise), ]
    list(d = d, n = n_testable)
  })
  # the CSV download resolves "<id>_df", so the table needs its own data-frame reactive
  sel_tbl_df <- reactive(sel_rank()$d)

  output$ptkey_tbl <- renderDT(datatable(
    unique(data.frame(`Manuscript patient` = PTKEY$patient,
                      `St Jude accession` = PTKEY$sj,
                      `Drug-response name` = PTKEY$bulk,
                      `Phylogeny / induction name` =
                        ifelse(PTKEY$id %in% c("4295","417","445","4084"), PTKEY$id, "-"),
                      `Genetic subtype` = PTKEY$subtype,
                      check.names = FALSE)),
    rownames = FALSE, options = list(pageLength = 8, dom = "tip")))

  output$sel_title <- renderText(sprintf(
    "Mutations ranked by rise under %s%s",
    if (identical(input$sel_drug, "DNR-Hi")) "daunorubicin" else "prednisolone",
    if (identical(input$sel_scope, "one")) sprintf(" - patient %s", input$dpat) else " - all patients"))

  output$sel_plot_title <- renderText(sprintf(
    "What did %s select for?%s",
    if (identical(input$sel_drug, "DNR-Hi")) "daunorubicin" else "prednisolone",
    if (identical(input$sel_scope, "one")) sprintf("  Patient %s.", input$dpat) else "  All five patients."))

  # One row per mutation: where it sat in DMSO, where it sat under the drug, and the
  # move between them. Reading a 92-row table to find that was the wrong ask.
  output$sel_plot <- renderPlot({
    d <- sel_tbl_df()
    validate(need(nrow(d) > 0, "No mutation was called under this drug in this selection."))
    d <- head(d[order(-d$Rise), ], 20)
    # gene AND the amino-acid change, or two variants in one gene collide
    d$lab <- make.unique(sprintf("%s  (pt %s)", d$Mutation, d$Patient))
    d$lab <- factor(d$lab, levels = rev(d$lab))
    d$evid <- ifelse(d$`Every replicate above`,
                     "every treated replicate above every control",
                     "replicates overlap the control")
    long <- rbind(
      data.frame(lab = d$lab, af = d$`DMSO %`,    what = "DMSO control", stringsAsFactors = FALSE),
      data.frame(lab = d$lab, af = d$`Treated %`, what = "Drug",         stringsAsFactors = FALSE))
    ggplot() +
      geom_segment(data = d, aes(y = lab, yend = lab, x = `DMSO %`, xend = `Treated %`,
                                 colour = evid),
                   arrow = arrow(length = unit(.16, "cm"), type = "closed"), linewidth = 1) +
      geom_point(data = long, aes(y = lab, x = af, shape = what), size = 2.6, colour = "#17212B") +
      scale_shape_manual(values = c(`DMSO control` = 1, Drug = 16), name = NULL) +
      scale_colour_manual(values = c(`every treated replicate above every control` = "#B03A2E",
                                     `replicates overlap the control` = "#B9C0C7"), name = NULL) +
      labs(x = "mutant allele frequency (%), DMSO to drug", y = NULL) +
      theme_lab(0) + theme(legend.position = "top", legend.box = "vertical")
  })

  output$sel_tbl <- renderDT(
    datatable(sel_tbl_df(), rownames = FALSE, options = list(pageLength = 10, dom = "ftip")))

  output$sel_note <- renderUI({
    s <- sel_rank(); n <- s$n
    hit <- if (nrow(s$d)) sum(s$d$`Every replicate above`) else 0
    note("Rise is mean treated minus mean DMSO for that patient, in allele-frequency ",
         "percent. With three treated against three control replicates, ", tags$b("every "),
         "treated replicate landing above every control replicate happens by chance with ",
         "p = 0.05, so across the ", n, " testable mutation-by-patient comparisons here ",
         "about ", round(n * 0.05), " would pass that filter with no drug effect at all; ",
         hit, " do. Read a single hit as a candidate to check, not as evidence of selection. ",
         "Variants flagged multi-mapped sit in paralogous loci where allele frequency is ",
         "unreliable, and synonymous changes are unlikely resistance drivers.")
  })

  output$dmut_rank_note <- renderUI({
    r <- sel_rank()$d
    if (!nrow(r)) return(note("No ranking for this drug."))
    i <- which(r$Mutation == input$dmut)
    drg <- if (identical(input$sel_drug, "DNR-Hi")) "daunorubicin" else "prednisolone"
    if (!length(i))
      return(note(sprintf("%s was not called in any patient under %s, so it is not in that ranking.",
                          input$dmut, drg)))
    row <- r[i[1], ]
    note(sprintf("Under %s in patient %s, %s goes from %.1f%% in DMSO to %.1f%% treated, a rise of %.1f points, ranked %d of %d. ",
                 drg, row$Patient, input$dmut, row$`DMSO %`, row$`Treated %`, row$Rise, i[1], nrow(r)),
         if (isTRUE(row$`Every replicate above`))
           "Every treated replicate exceeded every DMSO replicate - which about one comparison in twenty does by chance."
         else "Not every treated replicate exceeded every DMSO replicate.")
  })

  output$sj_plot <- renderPlot({
    d <- sj
    top <- names(sort(tapply(d$af, d$Mutation, max, na.rm = TRUE), decreasing = TRUE))[1:25]
    d <- d[d$Mutation %in% top, ]
    # Only 117 of the 25 x 9 combinations exist in the data. Left as-is, ggplot drops
    # them and the panel background shows through, which looks identical to a low
    # value. Completing the grid makes the gap explicit and colours it as missing.
    grid <- expand.grid(Treatment = sort(unique(d$Treatment)), Mutation = top,
                        stringsAsFactors = FALSE)
    d <- merge(grid, d[, c("Treatment", "Mutation", "af")], all.x = TRUE)
    d$Mutation <- factor(d$Mutation, levels = rev(top))
    ggplot(d, aes(Treatment, Mutation, fill = af)) +
      geom_tile(colour = "white", linewidth = .25) +
      scale_fill_viridis_c(option = "mako", direction = -1, name = "AF %",
                           na.value = "#E8ECEF") +
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
      scale_shape_manual(values = c(`Before induction` = 1, `After induction` = 16), name = NULL) +
      labs(x = "CD19", y = "CD34") + theme_lab()
  })

  # Raw counts would mislead: 30 cells were sequenced before induction against 85
  # after, so every clone looks larger afterwards. Each bar is the share of the
  # cells sequenced at that timepoint.
  clone_fig <- reactive({
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
  output$clone_plot <- renderPlot(print(clone_fig()))

  output$geno_plot <- renderPlot({
    ord <- order(cells$timepoint, ifelse(is.na(cells$clone), "zz", cells$clone))
    ids <- cells$Index[ord]
    ids <- ids[ids %in% rownames(GMAT)]
    m <- GMAT[ids, , drop = FALSE]
    d <- data.frame(cell = factor(rep(ids, ncol(m)), levels = ids),
                    variant = factor(rep(colnames(m), each = nrow(m)), levels = colnames(m)),
                    # states 1/2/4 are all mutant genotypes in the paper's own
                    # heatmap; only 0 is wild type
                    present = as.vector(m) > 0)
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
    sprintf("Paired before- and after-induction samples \u2014 patient 4295 only%s",
            if (identical(input$cpat, "4295")) "" else
              sprintf(" (you have patient %s selected above)", input$cpat)))

  output$qc_plot <- renderPlot({
    # n_snv drives the point size, so a row without it is dropped silently and the
    # panel comes out empty; require all three and say so when they are not there
    d <- cm(); d <- d[is.finite(d$ado) & is.finite(d$depth) & is.finite(d$n_snv), ]
    validate(need(nrow(d) > 0, sprintf(paste(
      "Depth and allelic dropout were not reported together with a mutation count",
      "for any single cell of patient %s, so there is nothing to plot here."),
      input$cpat)))
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

  # ---- carrier frequency before and after ----------------------------------
  pos_d <- reactive({
    d <- POS
    d$mutation <- d$gene
    dup <- d$gene[duplicated(d$gene)]
    d$mutation[d$gene %in% dup] <- paste(d$gene[d$gene %in% dup],
                                         sub(" .*$", "", d$locus[d$gene %in% dup]))
    d$fate <- ifelse(d$emergent, "Absent before induction",
              ifelse(d$delta > 0, "Rose", ifelse(d$delta < 0, "Fell", "Unchanged")))
    d[order(-d$pct_post), ]
  })

  output$pos_plot <- renderPlot({
    d <- pos_d()
    validate(need(nrow(d) > 0, "No genotype calls available."))
    long <- rbind(
      data.frame(mutation = d$mutation, fate = d$fate, when = "Before induction", pct = d$pct_pre),
      data.frame(mutation = d$mutation, fate = d$fate, when = "After induction",  pct = d$pct_post))
    long$when <- factor(long$when, levels = c("Before induction", "After induction"))
    lab <- d[order(-d$pct_post), ]
    lab$shifted <- stack_labels(lab$pct_post, 0, max(c(d$pct_pre, d$pct_post)), frac = .040)
    ggplot(long, aes(when, pct, group = mutation, colour = fate)) +
      geom_line(linewidth = .85, alpha = .85) +
      geom_point(size = 2.2) +
      geom_segment(data = lab, aes(x = 2, xend = 2.1, y = pct_post, yend = shifted),
                   linewidth = .3, colour = "#B9C0C7", show.legend = FALSE) +
      geom_text(data = lab, aes(x = 2.12, y = shifted, label = mutation),
                hjust = 0, size = 3.05, show.legend = FALSE) +
      scale_x_discrete(expand = expansion(mult = c(.08, .52))) +
      scale_y_continuous(labels = percent_format(accuracy = 1),
                         expand = expansion(mult = c(.03, .06))) +
      scale_colour_manual(values = c("Absent before induction" = CARD, "Rose" = TEAL,
                                     "Fell" = AMBER, "Unchanged" = GREY), name = NULL) +
      labs(x = NULL, y = "cells carrying the mutation") +
      theme_lab() + theme(legend.position = "top")
  })

  output$pos_head <- renderUI({
    d <- pos_d(); e <- d[d$emergent, ]
    note(sprintf(paste("%d variants across %d cells before induction and %d after.",
                       "%d were carried by no cell before induction and appear only after",
                       "induction: %s. %d rose, %d fell."),
                 nrow(d), d$cells_pre[1], d$cells_post[1], nrow(e),
                 if (nrow(e)) paste(sprintf("%s in %d of %d cells (%s)", e$gene, e$n_post,
                                            e$cells_post, percent(e$pct_post, accuracy = .1)),
                                    collapse = "; ") else "none",
                 sum(d$delta > 0 & !d$emergent), sum(d$delta < 0)))
  })

  pos_tbl_df <- reactive({
    d <- pos_d()
    data.frame(Gene = d$gene, Locus = d$locus,
                    `Cells before` = sprintf("%d / %d", d$n_pre, d$cells_pre),
                    `Cells after`  = sprintf("%d / %d", d$n_post, d$cells_post),
                    `% before` = round(100 * d$pct_pre, 1),
                    `% after`  = round(100 * d$pct_post, 1),
                    `Change (pp)` = round(100 * d$delta, 1),
               `Absent before` = ifelse(d$emergent, "yes", ""),
               check.names = FALSE)
  })
  output$pos_tbl <- renderDT(
    datatable(pos_tbl_df(), rownames = FALSE, options = list(pageLength = 10, dom = "ftip")))

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
      labs(x = "allele frequency before induction (sample 4272)",
           y = "allele frequency after induction (sample 4295)") + theme_lab()
  })

  af_tbl_df <- reactive({
    d <- AF[order(-AF$delta), c("gene", "locus", "pre", "post", "delta", "cells_pre", "cells_post")]
    names(d) <- c("Gene", "Locus", "Before", "After", "Change", "Cells before", "Cells after")
    d
  })
  output$af_tbl <- renderDT({
    datatable(af_tbl_df(), rownames = FALSE, options = list(pageLength = 8, dom = "ftip")) |>
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

  # ---- Induction therapy: the before/after exome phylogeny ----
  IND_TREE <- local({ t <- ape::read.tree(text = TREES[["clone"]]); t })
  IND_KEYS <- gsub("[._]", "-", IND_TREE$tip.label)

  # Per-tip composition. Signature fits exist only for the after-induction cells,
  # so the default is the cell's own alternate reads split by gene, which covers
  # 110 of the 113 tips including 29 of the 30 before-induction cells.
  ind_pie <- reactive({
    src <- SIG
    m <- matrix(0, nrow = length(IND_KEYS), ncol = ncol(src),
                dimnames = list(IND_KEYS, colnames(src)))
    i <- match(IND_KEYS, rownames(src)); ok <- !is.na(i)
    m[ok, ] <- as.matrix(src)[i[ok], , drop = FALSE]
    m[!is.finite(m)] <- 0
    tot <- colSums(m)
    keep <- names(sort(tot[tot > 0], decreasing = TRUE))
    keep <- keep[seq_len(min(8L, length(keep)))]
    rest <- setdiff(colnames(m)[colSums(m) > 0], keep)
    out <- m[, keep, drop = FALSE]
    if (length(rest)) out <- cbind(out, other = rowSums(m[, rest, drop = FALSE]))
    rs <- rowSums(out)
    out[rs > 0, ] <- out[rs > 0, , drop = FALSE] / rs[rs > 0]
    list(m = out, has = rs > 0)
  })

  ind_pal <- reactive({
    if (identical(input$ind_mark, "tp")) TP_EDGE
    else { pm <- ind_pie()$m
           setNames(grDevices::hcl.colors(ncol(pm), "Spectral"), colnames(pm)) }
  })

  # Clade membership straight off the topology. Node numbers are ape's internal
  # indices, which for this 113-tip tree run 114-135 and are the same numbers
  # Figure 7A prints above its pies.
  clade_tab <- reactive({
    t <- IND_TREE; nt <- ape::Ntip(t)
    tp <- as.character(TIPS$timepoint[match(IND_KEYS, TIPS$cell)])
    rows <- lapply(seq_len(t$Nnode) + nt, function(nd) {
      tips <- integer(0); stack <- nd
      while (length(stack)) {
        cur <- stack[1]; stack <- stack[-1]
        kids <- t$edge[t$edge[, 1] == cur, 2]
        tips <- c(tips, kids[kids <= nt]); stack <- c(stack, kids[kids > nt])
      }
      v <- tp[tips]
      data.frame(node = nd, n = length(tips),
                 pre  = sum(v == "Before induction", na.rm = TRUE),
                 post = sum(v == "After induction",  na.rm = TRUE))
    })
    d <- do.call(rbind, rows)
    d[d$n < nt, ]                      # drop the root, which is every cell
  })

  # ---- induction therapy, all four patients -------------------------------
  # The trees themselves are on the Phylogeny tab; this tab reads what induction
  # did to the mutations that sit on their branches.

  wgs_bulk <- reactive({
    q <- input$wgs_pt; req(q %in% WGS_PT)
    b <- BULK[BULK$patient == q, ]
    b <- b[is.finite(b$before) | is.finite(b$after), ]
    b$before[!is.finite(b$before)] <- 0; b$after[!is.finite(b$after)] <- 0
    b$mutation <- ifelse(nzchar(b$change), paste(b$gene, b$change), b$gene)
    b$fate <- ifelse(b$after > 0, "Still detectable", "Not detected after")
    b[order(-b$before), ]
  })


  # No panel should sit unchanged while a selector moves - a stale panel is
  # indistinguishable from a broken one. The single-cell exome experiment exists
  # for patient 4295 alone, so for any other patient these panels say so rather
  # than keep showing 4295.
  EXOME_PT <- "4295"
  exome_note <- function(q, what) div(
    class = "border rounded p-3 my-2", style = "background:#F7F9FA",
    tags$b(sprintf("Patient %s has no %s.", q, what)),
    div(style = "color:#5A6773;margin-top:.3rem",
        "The single-cell exome experiment, which sequenced cells from before and after ",
        "induction, was run on patient 4295 only. Select patient 4295 to see it."))
  exome_slot <- function(q, what, id, h)
    if (identical(as.character(q), EXOME_PT)) plotOutput(id, height = h) else exome_note(q, what)
  exome_tbl_slot <- function(q, what, id)
    if (identical(as.character(q), EXOME_PT)) tagList(DTOutput(id), dl_link(id)) else exome_note(q, what)


  ex <- function(x) identical(as.character(input$wgs_pt), EXOME_PT)
  output$ind_pies_head_slot <- renderUI(if (ex()) uiOutput("ind_pies_head"))
  output$pos_head_slot      <- renderUI(if (ex()) uiOutput("pos_head"))
  output$af_head_slot       <- renderUI(if (ex()) uiOutput("af_head"))

  # the descriptive text belongs to the 4295 phylogeny, so it goes with it
  output$ind_desc_slot <- renderUI({
    if (!identical(as.character(input$wgs_pt), EXOME_PT)) return(NULL)
    note("The Figure 7A single-cell phylogeny of patient 4295, spanning both samples: ",
         "4272 taken before induction and 4295 after four weeks of it. The branch into each ",
         "cell carries the sample it came from, and the pie at each internal node is the ",
         "before/after split of the cells beneath it, exactly as in Figure 7A. A node that ",
         "is entirely pink is a clade found only after induction. The pies are computed from ",
         "the topology here, rather than fixed to one tree as the figure script's hard-coded ",
         "clade vector is. The tree carries no branch lengths, so it is a cladogram: the ",
         "topology is meaningful and the horizontal distances are not. Which variants each ",
         "cell carries is the heatmap below, not the tip marks.")
  })
  output$clade_hdr  <- renderText(sprintf("Patient %s - before and after, clade by clade", input$wgs_pt))
  output$pospct_hdr <- renderText(sprintf("Patient %s - before and after induction, by percent of cells carrying the mutation", input$wgs_pt))
  output$postbl_hdr <- renderText(sprintf("Patient %s - every variant, by percent of cells", input$wgs_pt))
  output$afplot_hdr <- renderText(sprintf("Patient %s - allele frequency before and after induction", input$wgs_pt))
  output$aftbl_hdr  <- renderText(sprintf("Patient %s - every variant, before and after", input$wgs_pt))

  output$ind_legend_slot <- renderUI(
    if (identical(as.character(input$wgs_pt), EXOME_PT)) uiOutput("ind_legend"))
  output$ind_head_slot <- renderUI(
    if (identical(as.character(input$wgs_pt), EXOME_PT)) uiOutput("ind_head"))

  output$cell_plot_title <- renderText(sprintf(
    "Patient %s - single-cell genomes, before and after induction", input$cpat))
  output$clone_plot_title <- renderText(sprintf("Patient %s - clone composition", input$cpat))
  output$cell_plot_slot  <- renderUI(exome_slot(input$cpat, "paired single-cell exome samples", "cell_plot", 430))
  output$clone_plot_slot <- renderUI(exome_slot(input$cpat, "clone composition from paired exomes", "clone_plot", 430))

  output$ind_tree_slot   <- renderUI(exome_slot(input$wgs_pt, "single-cell exome phylogeny", "ind_tree", "auto"))
  output$ind_pies_slot   <- renderUI(exome_slot(input$wgs_pt, "clade composition", "ind_pies", 460))
  output$pos_plot_slot   <- renderUI(exome_slot(input$wgs_pt, "per-cell carrier frequencies", "pos_plot", 520))
  output$af_plot_slot    <- renderUI(exome_slot(input$wgs_pt, "single-cell allele frequencies", "af_plot", 470))
  output$geno_plot_slot  <- renderUI(exome_slot(input$wgs_pt, "single-cell genotype matrix", "geno_plot", 420))
  output$clone_plot2_slot<- renderUI(exome_slot(input$wgs_pt, "clone composition", "clone_plot2", 430))
  output$pos_tbl_slot    <- renderUI(exome_tbl_slot(input$wgs_pt, "per-cell carrier frequencies", "pos_tbl"))
  output$af_tbl_slot     <- renderUI(exome_tbl_slot(input$wgs_pt, "single-cell allele frequencies", "af_tbl"))
  output$geno_title  <- renderText(sprintf("Patient %s - genotypes across 31 variants", input$wgs_pt))
  output$clone2_title<- renderText(sprintf("Patient %s - clone composition before and after", input$wgs_pt))

  output$wgs_card_title <- renderText(
    sprintf("Patient %s - before and after induction, measured in the bulk", input$wgs_pt))
  output$wgs_tbl_title <- renderText(
    sprintf("Patient %s - every branch mutation, before and after", input$wgs_pt))

  output$wgs_slope <- renderPlot({
    b <- wgs_bulk()
    validate(need(nrow(b) > 0, "No bulk measurement for this patient."))
    long <- rbind(
      data.frame(mutation = b$mutation, fate = b$fate, when = "Before induction", vaf = b$before),
      data.frame(mutation = b$mutation, fate = b$fate, when = "After induction",  vaf = b$after))
    long$when <- factor(long$when, levels = c("Before induction", "After induction"))
    # ggrepel is not in the webR repo, so labels are de-overlapped here: work in
    # the sqrt space the axis actually uses, then push each label down until it
    # clears the one above by a fixed fraction of the panel.
    lab <- b[b$after > 0, ]
    lab <- lab[order(-lab$after), ]
    if (nrow(lab)) {
      # the axis is sqrt-transformed, so space the labels in sqrt units
      top <- sqrt(max(c(b$before, b$after), na.rm = TRUE))
      lab$shifted <- stack_labels(sqrt(lab$after), 0, top, frac = .042)^2
    }
    ggplot(long, aes(when, vaf, group = mutation, colour = fate)) +
      geom_line(linewidth = .8, alpha = .8) +
      geom_point(size = 2.1) +
      geom_segment(data = lab, aes(x = 2, xend = 2.12, y = after, yend = shifted),
                   linewidth = .3, colour = "#B9C0C7", show.legend = FALSE) +
      geom_text(data = lab, aes(x = 2.14, y = shifted, label = mutation),
                hjust = 0, size = 3.05, show.legend = FALSE) +
      scale_x_discrete(expand = expansion(mult = c(.08, .62))) +
      scale_y_sqrt(labels = percent_format(accuracy = 1),
                   breaks = c(0, .01, .05, .1, .2, .3, .4, .5, .6)) +
      scale_colour_manual(values = c("Still detectable" = CARD,
                                     "Not detected after" = "#9AA3AB"), name = NULL) +
      labs(x = NULL, y = "Variant allele frequency in bulk") +
      theme_lab() + theme(legend.position = "top")
  })

  output$wgs_slope_head <- renderUI({
    b <- wgs_bulk()
    kept <- b[b$after > 0, ]
    note(sprintf(paste("%d branch mutations. Median allele frequency fell from %s at diagnosis",
                       "to %s in remission. %d were still detectable in the remission bulk%s."),
                 nrow(b), percent(median(b$before), accuracy = .1),
                 percent(median(b$after), accuracy = .1), nrow(kept),
                 if (nrow(kept)) sprintf(", the highest being %s at %s",
                                         kept$mutation[which.max(kept$after)],
                                         percent(max(kept$after), accuracy = .1)) else ""))
  })

  wgs_tbl_df <- reactive({
    b <- wgs_bulk()
    data.frame(Mutation = b$mutation, Effect = b$effect, Locus = b$locus,
                    `Before induction` = round(b$before, 4), `After induction` = round(b$after, 4),
                    `Fold change` = ifelse(b$after > 0, round(b$before / b$after, 1), NA),
               `Depth before` = b$dp_before, `Depth after` = b$dp_after,
               check.names = FALSE)
  })
  output$wgs_tbl <- renderDT(
    datatable(wgs_tbl_df(), rownames = FALSE, options = list(pageLength = 10, dom = "ftip")))

  output$ind_pies <- renderPlot({
    d <- clade_tab()
    d <- d[d$n >= input$ind_min & (d$pre + d$post) > 0, ]
    validate(need(nrow(d) > 0, "No clade reaches that size."))
    d <- d[order(-d$n), ]
    lab <- setNames(sprintf("node %d  (%d cells)", d$node, d$n), d$node)
    long <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) data.frame(
      node = factor(unname(lab[as.character(d$node[i])]), levels = unname(lab)),
      sample = factor(c("Before induction", "After induction"),
                      levels = c("Before induction", "After induction")),
      frac = c(d$pre[i], d$post[i]) / (d$pre[i] + d$post[i]))))
    ggplot(long, aes(x = "", y = frac, fill = sample)) +
      geom_col(width = 1, colour = "white", linewidth = .5) +
      coord_polar(theta = "y", direction = -1) +
      facet_wrap(~ node, ncol = 6) +
      scale_fill_manual(values = TP_EDGE, name = NULL) +
      theme_void(base_size = 12) +
      theme(legend.position = "top",
            strip.text = element_text(size = 10, margin = margin(2, 2, 4, 2)))
  })

  output$ind_pies_head <- renderUI({
    d <- clade_tab(); d <- d[d$n >= input$ind_min, ]
    only_post <- d[d$pre == 0 & d$post > 0, ]
    only_pre  <- d[d$post == 0 & d$pre > 0, ]
    note(sprintf(paste("%d clades of at least %d cells. %d hold only cells from after induction",
                       "(%s) and %d only cells from before it%s."),
                 nrow(d), input$ind_min, nrow(only_post),
                 if (nrow(only_post)) paste("nodes", paste(only_post$node, collapse = ", ")) else "none",
                 nrow(only_pre),
                 if (nrow(only_pre)) sprintf(" (nodes %s)", paste(only_pre$node, collapse = ", ")) else ""))
  })

  output$ind_title <- renderText(
    if (!identical(as.character(input$wgs_pt), EXOME_PT))
      sprintf("Single-cell exome phylogeny \u2014 not available for patient %s", input$wgs_pt) else
    sprintf("Patient 4295, before and after induction \u2014 %d single-cell genomes",
            ape::Ntip(IND_TREE)))

  output$ind_tree <- renderPlot(height = function() {
    per <- if (identical(input$ind_mark, "tp")) 11 else 30
    max(560, min(3600, round(per * ape::Ntip(IND_TREE))))
  }, {
    t <- IND_TREE
    tp <- TIPS$timepoint[match(IND_KEYS, TIPS$cell)]
    ecol <- rep("#5A6570", nrow(t$edge)); ewid <- rep(.9, nrow(t$edge))
    term <- t$edge[, 2] <= ape::Ntip(t)
    v <- as.character(tp)[t$edge[term, 2]]
    ecol[term] <- ifelse(is.na(v), "#C7CDD2", unname(TP_EDGE[v]))
    ewid[term] <- 2.1
    # plot.phylo sizes the x axis to the tree alone, so on a cladogram every tip
    # sits exactly on the right edge and half of each signature disc is cut off.
    # Measure the axis first, then re-plot with headroom for the discs.
    par(mar = c(1, 1, 1, 1), xpd = TRUE)
    xl <- NULL
    if (!identical(input$ind_mark, "tp"))
      xl <- tryCatch({
        pp <- plot(t, type = input$ind_type, show.tip.label = isTRUE(input$ind_lab),
                   cex = .6, use.edge.length = FALSE, plot = FALSE)
        c(pp$x.lim[1], pp$x.lim[1] + diff(pp$x.lim) * 1.10)
      }, error = function(e) NULL)
    plot(t, type = input$ind_type, show.tip.label = isTRUE(input$ind_lab),
         cex = .6, no.margin = FALSE, edge.color = ecol, edge.width = ewid,
         use.edge.length = FALSE, x.lim = xl)
    # Figure 7A draws a pie at each major clade giving the pre/post split of the
    # cells beneath it. Computed here from the topology rather than hard-coded to
    # one tree, so it survives the tree being rebuilt.
    if (isTRUE(input$ind_node) && identical(input$ind_type, "phylogram")) {
      nt <- ape::Ntip(t)
      desc <- lapply(seq_len(t$Nnode) + nt, function(nd) {
        tips <- integer(0); stack <- nd
        while (length(stack)) {
          cur <- stack[1]; stack <- stack[-1]
          kids <- t$edge[t$edge[, 1] == cur, 2]
          tips <- c(tips, kids[kids <= nt]); stack <- c(stack, kids[kids > nt])
        }
        tips
      })
      sizes <- lengths(desc)
      sel <- which(sizes >= input$ind_min & sizes < nt)
      if (length(sel)) {
        pm <- t(vapply(desc[sel], function(ti) {
          v <- as.character(tp)[ti]
          c(`Before induction` = sum(v == "Before induction", na.rm = TRUE),
            `After induction`  = sum(v == "After induction",  na.rm = TRUE))
        }, numeric(2)))
        rs <- rowSums(pm); ok <- rs > 0
        pm[ok, ] <- pm[ok, , drop = FALSE] / rs[ok]
        if (any(ok))
          ape::nodelabels(pie = pm[ok, , drop = FALSE], node = (sel + nt)[ok],
                          piecol = unname(TP_EDGE), cex = .45)
      }
    }
    if (identical(input$ind_mark, "tp")) {
      cols <- ifelse(is.na(tp), GREY, unname(TP_EDGE[as.character(tp)]))
      ape::tiplabels(pch = 19, col = cols, cex = 1.1)
    } else {
      pq <- ind_pie(); pm <- pq$m
      per_tip <- (par("din")[2] * 96) / max(ape::Ntip(t), 1)
      pc <- 0.70 * min(1, per_tip / 30)
      if (any(!pq$has)) ape::tiplabels(pch = 19, col = GREY, cex = .6, tip = which(!pq$has))
      if (any(pq$has))  ape::tiplabels(pie = pm[pq$has, , drop = FALSE], tip = which(pq$has),
                                       piecol = unname(ind_pal()), cex = pc)
    }
  })

  output$ind_legend <- renderUI({
    dot <- function(c) span(style = sprintf("display:inline-block;width:.75rem;height:.75rem;border-radius:50%%;background:%s;margin-right:.3rem;vertical-align:-1px", c))
    bar <- function(c) span(style = sprintf("display:inline-block;width:1.1rem;height:.2rem;background:%s;margin-right:.3rem;vertical-align:.18rem", c))
    pal <- ind_pal()
    tagList(
      div(style = "display:flex;flex-wrap:wrap;gap:.35rem 1rem;margin:.4rem 0 0;font-size:.85rem",
          span(style = "color:#64707C",
               if (identical(input$ind_mark, "tp")) "tips:"
               else if (identical(input$ind_mark, "sig")) "signature:" else "gene:"),
          lapply(names(pal), function(n) span(dot(pal[[n]]), n))),
      div(style = "display:flex;flex-wrap:wrap;gap:.35rem 1rem;margin:.3rem 0 0;font-size:.85rem",
          span(style = "color:#64707C", "terminal branch and clade pie:"),
          lapply(names(TP_EDGE), function(n) span(bar(TP_EDGE[[n]]), n))))
  })

  output$ind_head <- renderUI({
    pq <- ind_pie()
    tp <- TIPS$timepoint[match(IND_KEYS, TIPS$cell)]
    miss <- sum(!pq$has)
    bits <- sprintf("%d cells: %d before induction, %d after",
                    length(IND_KEYS), sum(tp == "Before induction", na.rm = TRUE),
                    sum(tp == "After induction", na.rm = TRUE))
    if (!identical(input$ind_mark, "tp"))
      bits <- c(bits, sprintf("%d %s no %s and are drawn as a plain dot", miss,
                              if (miss == 1) "tip has" else "tips have",
                              "signature fit (no cell from before induction was fitted)"))
    bits <- c(bits, sprintf("clones %s appear only after induction",
                            paste(EMERGENT, collapse = " and ")))
    note(paste0(paste(bits, collapse = "; "), "."))
  })

  output$clone_plot2 <- renderPlot(print(clone_fig()))

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
    # the per-cell signature pie is the default wherever the tree has signature
    # fits; pick it by name, since which option sits at o[[2]] depends on which
    # measurements that particular tree happens to carry
    sel <- if ("pie" %in% o) "pie" else if (length(o) > 1) o[[2]] else o[[1]]
    radioButtons("tipcol", "Colour tips by", choices = o, selected = sel)
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

  # the figure's own palette: 4272 = pre, 4295 = post
  TP_EDGE <- c(`Before induction` = "#56B4E9", `After induction` = "#CC79A7")

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
        c(`Before induction` = "#4A555F", `After induction` = CARD, `not annotated` = GREY)[lv]
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

  # Exome-mapped mutations placed on the branches of the tree being drawn. Only
  # the four Figure 5 patient trees carry them; the in vitro benchmark and the
  # Figure 7 clone tree get nothing.
  tree_ann <- reactive({
    q <- input$tree
    empty <- data.frame(node = integer(0), n = integer(0), label = character(0),
                        kept = logical(0), tag = character(0))
    if (!isTRUE(input$treemut) || is.null(WGSMUT[[q]])) return(empty)
    m <- WGSMUT[[q]]$mut
    if (!nrow(m)) return(empty)
    t <- cur_tree(); all <- sort(t$tip.label)
    ck <- vapply(node_sets(t), split_key, character(1), all = all)
    ck[!(seq_along(ck) %in% t$edge[, 2])] <- NA_character_   # the root has no branch
    m$node <- match(vapply(strsplit(m$key, "|", fixed = TRUE), split_key,
                           character(1), all = all), ck)
    m <- m[!is.na(m$node), ]
    if (!nrow(m)) return(empty)
    b <- BULK[BULK$patient == q, ]
    genes <- lapply(strsplit(m$label, "; *"), function(x) sub(" .*$", "", trimws(x)))
    m$kept <- vapply(genes, function(g) {
      v <- b$after[b$gene %in% g]
      length(v) > 0 && any(is.finite(v) & v > 0)
    }, logical(1))
    m <- m[order(-m$n), ]
    m$tag <- paste0("M", seq_len(nrow(m)))
    m
  })

  # Sixteen gene names cannot be written along a branch without burying the tree,
  # so the branch carries a short tag and the tags are expanded underneath.
  output$tree_muts <- renderUI({
    a <- tree_ann()
    if (!nrow(a)) return(NULL)
    rows <- lapply(seq_len(nrow(a)), function(i) tags$li(
      style = "margin-bottom:.25rem",
      tags$b(style = sprintf("color:%s", if (a$kept[i]) CARD else "#39424A"), a$tag[i]),
      sprintf(" \u00b7 %s \u00b7 ", if (a$n[i] == 1) "1 cell" else sprintf("%d cells", a$n[i])),
      a$label[i],
      tags$span(class = "text-muted",
                if (a$kept[i]) " \u2014 still detectable in the remission bulk"
                else " \u2014 not detected in the remission bulk")))
    div(style = "font-size:.86rem;margin-top:.6rem",
        div(class = "text-muted", style = "margin-bottom:.3rem",
            "Mutations mapped onto branches. A tag on the trunk is carried by every ",
            "sampled cell; a tag further out belongs to one clone or one cell. Red means ",
            "the mutation was still detectable in that patient's remission bulk sample."),
        tags$ul(style = "padding-left:1.1rem", rows))
  })

  output$tree_plot <- renderPlot(height = function() {
    n <- tryCatch(ape::Ntip(cur_tree()), error = function(e) 40)
    per <- if (identical(input$tipcol, "pie")) 30 else 9   # pies need room to read
    max(560, min(3600, round(per * n)))
  }, {
    t <- cur_tree(); ann <- tip_ann()
    lv <- sort(unique(ann))
    pal <- if (identical(input$tipcol, "timepoint"))
             c(`Before induction` = "#4A555F", `After induction` = CARD, `not annotated` = GREY)[lv]
           else if (input$tipcol %in% c("ado", "depth")) {
             o <- lv[order(suppressWarnings(as.numeric(sub("^[\\[(]", "", sub(",.*", "", lv)))))]
             setNames(grDevices::hcl.colors(length(o), "Viridis"), o)[lv]
           } else setNames(grDevices::hcl.colors(length(lv), "Spectral"), lv)
    pal[is.na(pal)] <- GREY
    cols <- pal[ann]
    # The terminal branch of each cell is coloured by the sample that cell came
    # from, so before and after induction are readable straight off the topology.
    tp <- tip_timepoint()
    ecol <- rep("#5A6570", nrow(t$edge)); ewid <- rep(.9, nrow(t$edge))
    if (!is.null(tp)) {
      term <- t$edge[, 2] <= ape::Ntip(t)
      v <- tp[t$edge[term, 2]]
      ecol[term] <- ifelse(is.na(v), "#C7CDD2", unname(TP_EDGE[v]))
      ewid[term] <- 2.1
    }
    ann_m <- tree_ann()
    if (nrow(ann_m)) {
      e <- match(ann_m$node, t$edge[, 2]); ok <- !is.na(e)
      ecol[e[ok]] <- ifelse(ann_m$kept[ok], CARD, "#39424A")
      ewid[e[ok]] <- 2.6
    }
    par(mar = c(1, 1, 1, 1), xpd = TRUE)
    plot(t, type = input$ttype, show.tip.label = isTRUE(input$tiplab),
         tip.color = cols, cex = .6, no.margin = FALSE,
         edge.color = ecol, edge.width = ewid,
         use.edge.length = !is.null(t$edge.length))
    if (nrow(ann_m)) {
      e <- match(ann_m$node, t$edge[, 2]); ok <- !is.na(e)
      ape::edgelabels(ann_m$tag[ok], e[ok], frame = "rect", cex = .68, adj = c(.5, .5),
                      bg = ifelse(ann_m$kept[ok], "#F4DADA", "#EDF0F2"),
                      col = ifelse(ann_m$kept[ok], CARD, "#39424A"), font = 2)
    }
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

  gene_tbl_df <- reactive({
    d <- gene_am()[, c("gene", "label", "nonsyn", "fold_size_corr", "relapse_freq_pct", "mean_am")]
    names(d) <- c("Gene", "Class", "Nonsyn. events", "Fold recurrence", "Relapse %", "Mean AlphaMissense")
    d[order(-d$`Fold recurrence`), ]
  })
  output$gene_tbl <- renderDT({
    datatable(gene_tbl_df(), rownames = FALSE, selection = "single",
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
