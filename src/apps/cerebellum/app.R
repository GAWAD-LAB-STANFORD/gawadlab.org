# Developing Cerebellum Atlas browser
# ---------------------------------------------------------------------------
# Carter RA et al., Current Biology 28, 2910-2920 (2018): single-cell RNA-seq of
# the developing mouse cerebellum, E10 to P10. Map any gene onto the embedding,
# follow it through development, and rank genes by cell-type specificity.
# Runs entirely in the browser - no server, no install.
#
#   Rscript src/apps/cerebellum/prep_bundle.R DATA_DIR bundle.rds
#   shiny::runApp("src/apps/cerebellum")
# ---------------------------------------------------------------------------

library(shiny); library(bslib); library(ggplot2); library(DT)

B <- local({
  for (p in c("bundle.rds", file.path("..", "bundle.rds"),
              file.path("src", "apps", "cerebellum", "bundle.rds")))
    if (file.exists(p)) return(readRDS(p))
  stop("bundle.rds not found - build it with prep_bundle.R first")
})

cells <- B$cells; bins <- B$bins; gsz <- B$gsz
GENES <- B$genes                      # genes with a position on the embedding
GMAT  <- B$gmat; SCALE <- B$vmax / 255
CT <- B$CT; TP <- B$TP
CT_MEAN <- B$CT_MEAN / B$MSCALE; CT_PCT <- B$CT_PCT / B$PSCALE
TP_MEAN <- B$TP_MEAN / B$MSCALE; TP_PCT <- B$TP_PCT / B$PSCALE
SUMM <- B$summ                        # every gene with group statistics
N_CT <- setNames(gsz$n_cells[gsz$kind == "cell_type"], gsz$group[gsz$kind == "cell_type"])

CT_COL <- c(Granule = "#2F5D70", Progenitor = "#F6A30C", Purkinje = "#8C1515",
            GABA = "#5B8C5A", Glia = "#7E57C2", NTZ = "#C96A2B",
            Interneuron = "#3E8FA8", Other = "#B9C0C7")[CT]
CT_COL[is.na(CT_COL)] <- "#B9C0C7"; names(CT_COL) <- CT
TP_COL <- setNames(grDevices::hcl.colors(length(TP), "Zissou 1"), TP)

# Tile size is the real spacing between neighbouring grid centres, with a hair of
# overlap: an undersized tile leaves a sub-pixel seam at every edge, and those
# seams alias into white lines ruled across the map.
step <- function(v) { u <- sort(unique(round(v, 6))); if (length(u) > 1) min(diff(u)) else 1 }
BINW <- step(bins$x) * 1.02; BINH <- step(bins$y) * 1.02

theme_lab <- function(rot = 0) {
  theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          axis.text.x = element_text(angle = rot, hjust = if (rot > 0) 1 else .5))
}
note <- function(...) div(class = "text-muted", style = "font-size:.85rem;margin-top:.5rem", ...)
ext  <- function(label, url) a(label, href = url, target = "_blank", rel = "noopener")
gene_links <- function(g) tagList(
  ext("CZ CELLxGENE", paste0("https://cellxgene.cziscience.com/gene-expression?genes=", g)), " · ",
  ext("MGI", paste0("https://www.informatics.jax.org/quicksearch/summary?query=", g)), " · ",
  ext("Allen Developing Mouse Brain", paste0("https://developingmouse.brain-map.org/search/show?search_term=", g)))
elsewhere <- function(g, what, why = NULL)
  div(class = "border rounded p-3", style = "background:#F7F8F9;font-size:.9rem",
      div(strong(g), sprintf(" has no %s in this atlas.", what)),
      if (!is.null(why)) div(class = "text-muted", style = "margin-top:.25rem", why),
      div(style = "margin-top:.6rem", "Look it up in ", gene_links(g), "."))

gvec <- function(g) { i <- match(g, GENES)
                      if (is.na(i)) NULL else as.integer(GMAT[, i]) * SCALE }
# Squares where the gene was never detected are drawn grey rather than at the
# bottom of the colour ramp, so "not detected" cannot be read as "low", and the
# ramp is square-root: most genes sit at zero in most squares, and a linear ramp
# turns a marker confined to one population into a flat wash.
emap <- function(g, title = NULL) {
  v <- gvec(g); req(!is.null(v))
  d <- bins; d$mean <- v; z <- d$mean <= 0
  ggplot() +
    geom_tile(data = d[z, ], aes(x, y), fill = "#EDF0F2", width = BINW, height = BINH) +
    geom_tile(data = d[!z, ], aes(x, y, fill = mean), width = BINW, height = BINH) +
    scale_fill_viridis_c(option = "rocket", direction = -1, trans = "sqrt",
                         name = sprintf("%s\nmean log expr", g)) +
    labs(x = "t-SNE 1", y = "t-SNE 2", title = title) +
    coord_equal() + theme_lab() + theme(panel.grid = element_blank())
}

SCOPES <- list("Every gene" = SUMM$gene, "Mapped on the atlas" = sort(GENES))
SCOPE_LABELS <- sprintf("%s (%s)", names(SCOPES),
                        format(lengths(SCOPES), big.mark = ",", trim = TRUE))

# --- ui --------------------------------------------------------------------
ui <- page_navbar(
  id = "nav",
  title = "Developing Cerebellum Atlas",
  theme = bs_theme(version = 5, primary = "#2F5D70"),
  header = tags$style(HTML(
    ".navbar .navbar-brand{font-size:1.45rem;font-weight:700}",
    ".navbar .nav-link{font-size:1.18rem;font-weight:600;padding:.5rem 1rem}",
    ".navbar .nav-link.active{font-weight:700}")),
  sidebar = sidebar(
    width = 300,
    selectizeInput("gene", "Selected gene", choices = SCOPES[[1]],
                   selected = if ("Atoh1" %in% SUMM$gene) "Atoh1" else SUMM$gene[1],
                   options = list(placeholder = "type a gene, e.g. Atoh1", maxOptions = 200)),
    radioButtons("glist", "Limit the list to", choiceNames = SCOPE_LABELS,
                 choiceValues = names(SCOPES), selected = names(SCOPES)[1]),
    conditionalPanel("input.nav == 'Compare'", hr(),
      selectizeInput("gene2", "Second gene", choices = sort(GENES),
                     selected = if ("Ptf1a" %in% GENES) "Ptf1a" else GENES[2],
                     options = list(placeholder = "type a gene, e.g. Ptf1a", maxOptions = 200))),
    uiOutput("gene_card")
  ),

  nav_panel("Atlas",
    card(card_header("Developing mouse cerebellum, E10 to P10"),
      layout_columns(col_widths = c(3, 9),
        radioButtons("fill", "Colour by",
                     c("Gene expression" = "expr", "Cell type" = "type",
                       "Timepoint" = "time", "Cell density" = "dens")),
        div(uiOutput("emb_miss"), plotOutput("emb", height = 540))),
      note(sprintf("%s cells on the published t-SNE, across %d annotated populations and %d collection timepoints. ",
                   format(nrow(cells), big.mark = ","), length(CT), length(TP)),
           "Gene expression is the mean log expression of every cell in each grid square, ",
           sprintf("so a sparsely populated square is not over-read. %s genes have enough cells to map.",
                   format(length(GENES), big.mark = ","))))),

  nav_panel("Gene profile",
    layout_columns(col_widths = c(6, 6),
      card(card_header(textOutput("ct_title")), plotOutput("ct_plot", height = 380),
           note("The mean runs over every cell in the type, so cells with no detected transcript ",
                "pull it down. Point size is the percentage of cells with any detected transcript.")),
      card(card_header("Expression through development"), plotOutput("tp_plot", height = 380),
           note("One point per collection timepoint, embryonic day 10 to postnatal day 10."))),
    card(card_header("Where it sits on the atlas"), uiOutput("prof_miss"),
         plotOutput("prof_map", height = 460))),

  nav_panel("Markers",
    card(card_header("Rank genes by cell-type specificity"),
      layout_columns(col_widths = c(3, 9),
        div(selectInput("mk_type", "Cell type", c("Any", CT)),
            sliderInput("mk_pct", "Min % of cells detected", 0, 100, 10, step = 5),
            sliderInput("mk_spec", "Min specificity", 0, 1, 0.3, step = 0.05),
            checkboxInput("mk_mapped", "Only genes mapped on the atlas", TRUE),
            downloadButton("mk_dl", "Download (CSV)", class = "btn-sm btn-primary")),
        DTOutput("mk_tbl")),
      note("Specificity is the peak cell type's mean expression as a share of the sum across all ",
           sprintf("%d types: 1.0 is seen in one type only, %.3f is flat across every type. ", length(CT), 1 / length(CT)),
           "Fold over 2nd compares the peak type with the runner-up. Click a row to select that gene."))),

  nav_panel("Compare",
    layout_columns(col_widths = c(6, 6),
      card(plotOutput("cmp_a", height = 420)), card(plotOutput("cmp_b", height = 420))),
    card(card_header("Per-region comparison"), uiOutput("cmp_miss"), plotOutput("cmp_sc", height = 400),
         note("Each point is one grid square, coloured by the cell type that dominates it and sized ",
              "by how many cells it holds. Squares off the diagonal are where the two genes disagree.")))
)

# --- server ----------------------------------------------------------------
server <- function(input, output, session) {
  g  <- reactive(if (is.null(input$gene)  || !nzchar(input$gene))  SUMM$gene[1] else input$gene)
  g2 <- reactive(if (is.null(input$gene2) || !nzchar(input$gene2)) GENES[1]     else input$gene2)
  mapped  <- reactive(g()  %in% GENES)
  mapped2 <- reactive(g2() %in% GENES)
  row_of  <- function(x) SUMM[match(x, SUMM$gene), ]

  observeEvent(input$glist, {
    ch <- SCOPES[[input$glist]]; cur <- isolate(input$gene)
    updateSelectizeInput(session, "gene", choices = ch,
                         selected = if (!is.null(cur) && cur %in% ch) cur else ch[1])
  }, ignoreInit = TRUE)

  NOMAP <- "Genes detected in too few cells to bin reliably are left off the map rather than drawn from noise."
  output$emb_miss  <- renderUI(
    if (identical(input$fill, "expr") && !mapped()) elsewhere(g(), "position on the embedding", NOMAP))
  output$prof_miss <- renderUI(if (!mapped())  elsewhere(g(),  "position on the embedding", NOMAP))
  output$cmp_miss  <- renderUI(
    if (!mapped() || !mapped2())
      elsewhere(if (!mapped()) g() else g2(), "position on the embedding", NOMAP))

  output$gene_card <- renderUI({
    r <- row_of(g()); req(!is.na(r$gene))
    i <- match(g(), rownames(CT_PCT))
    tagList(hr(),
      h4(r$gene, style = "margin:.2rem 0 0"),
      div(class = "text-muted", style = "font-size:.85rem",
          sprintf("peaks in %s at %s", r$peak_type, r$peak_time)),
      tags$table(class = "table table-sm", style = "font-size:.85rem;margin-top:.5rem",
        tags$tbody(
          tags$tr(tags$td("Peak cell type"), tags$td(r$peak_type)),
          tags$tr(tags$td("Specificity"),    tags$td(sprintf("%.2f", r$specificity))),
          tags$tr(tags$td("Fold over 2nd"),  tags$td(sprintf("%.1f×", r$fold_2nd))),
          tags$tr(tags$td("Peak timepoint"), tags$td(r$peak_time)),
          tags$tr(tags$td("Detected in"),    tags$td(sprintf("%.1f%% of %s cells",
                                                    r$max_pct, r$peak_type))))),
      if (!r$mapped) note("Not mapped on the atlas — too few cells for a reliable grid."),
      div(style = "font-size:.8rem;margin-top:.5rem", "Look up ", strong(r$gene), " in ", gene_links(r$gene)))
  })

  output$emb <- renderPlot({
    switch(input$fill,
      expr = { req(mapped()); emap(g()) },
      type = ggplot(cells, aes(D1, D2, colour = cell_type)) +
        geom_point(size = .35, alpha = .5) +
        scale_colour_manual(values = CT_COL, name = NULL) +
        guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
        labs(x = "t-SNE 1", y = "t-SNE 2") + coord_equal() + theme_lab(),
      time = ggplot(cells, aes(D1, D2, colour = timepoint)) +
        geom_point(size = .35, alpha = .5) +
        scale_colour_manual(values = TP_COL, name = NULL) +
        guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
        labs(x = "t-SNE 1", y = "t-SNE 2") + coord_equal() + theme_lab(),
      dens = ggplot(bins, aes(x, y, fill = n_cells)) +
        geom_tile(width = BINW, height = BINH) +
        scale_fill_viridis_c(option = "mako", direction = -1, trans = "sqrt", name = "cells") +
        labs(x = "t-SNE 1", y = "t-SNE 2") + coord_equal() + theme_lab() +
        theme(panel.grid = element_blank()))
  })

  output$ct_title <- renderText(sprintf("%s across cell types", g()))
  output$ct_plot <- renderPlot({
    i <- match(g(), rownames(CT_MEAN)); req(!is.na(i))
    d <- data.frame(type = factor(CT, levels = CT[order(CT_MEAN[i, ])]),
                    mean = CT_MEAN[i, ], pct = CT_PCT[i, ])
    ggplot(d, aes(mean, type)) +
      geom_segment(aes(x = 0, xend = mean, yend = type), colour = "#D6DBE0", linewidth = 1) +
      geom_point(aes(size = pct, colour = type)) +
      scale_colour_manual(values = CT_COL, guide = "none") +
      scale_size_continuous(range = c(2, 9), name = "% detected", limits = c(0, NA)) +
      labs(x = "mean log expression (all cells in the type)", y = NULL) + theme_lab()
  })

  output$tp_plot <- renderPlot({
    i <- match(g(), rownames(TP_MEAN)); req(!is.na(i))
    d <- data.frame(tp = factor(TP, levels = TP), mean = TP_MEAN[i, ], pct = TP_PCT[i, ])
    ggplot(d, aes(tp, mean, group = 1)) +
      geom_line(colour = "#2F5D70", linewidth = 1) +
      geom_point(aes(size = pct), colour = "#2F5D70") +
      scale_size_continuous(range = c(2, 8), name = "% detected", limits = c(0, NA)) +
      labs(x = NULL, y = "mean log expression", title = g()) + theme_lab(45)
  })

  output$prof_map <- renderPlot({ req(mapped()); emap(g()) })

  mk <- reactive({
    d <- SUMM
    if (isTRUE(input$mk_mapped)) d <- d[d$mapped, ]
    if (!identical(input$mk_type, "Any")) d <- d[d$peak_type == input$mk_type, ]
    d <- d[d$max_pct >= input$mk_pct & d$specificity >= input$mk_spec, ]
    d[order(-d$specificity, -d$max_pct), c("gene", "peak_type", "specificity",
                                           "fold_2nd", "peak_time", "max_pct")]
  })
  output$mk_tbl <- renderDT(
    datatable(mk(), rownames = FALSE, selection = "single",
              colnames = c("Gene", "Peak cell type", "Specificity", "Fold over 2nd",
                           "Peak timepoint", "Max % detected"),
              options = list(pageLength = 15, dom = "ftip", order = list())))
  observeEvent(input$mk_tbl_rows_selected, {
    i <- input$mk_tbl_rows_selected
    if (length(i)) updateSelectizeInput(session, "gene", selected = mk()$gene[i])
  })
  output$mk_dl <- downloadHandler(
    function() sprintf("cerebellum_markers_%s.csv", Sys.Date()),
    function(f) write.csv(mk(), f, row.names = FALSE))

  output$cmp_a <- renderPlot({ req(mapped());  emap(g(),  g())  })
  output$cmp_b <- renderPlot({ req(mapped2()); emap(g2(), g2()) })
  output$cmp_sc <- renderPlot({
    req(mapped(), mapped2())
    d <- data.frame(a = gvec(g()), b = gvec(g2()), type = bins$cell_type, n = bins$n_cells)
    ggplot(d, aes(a, b, colour = type, size = n)) +
      geom_abline(slope = 1, intercept = 0, colour = "#C7CDD2", linetype = 2) +
      geom_point(alpha = .7) +
      scale_colour_manual(values = CT_COL, name = NULL) +
      scale_size_continuous(range = c(.6, 5), guide = "none") +
      guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
      labs(x = sprintf("%s mean log expression", g()),
           y = sprintf("%s mean log expression", g2())) + theme_lab()
  })
}

shinyApp(ui, server)
