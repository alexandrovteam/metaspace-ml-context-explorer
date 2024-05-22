
# Libraries ---------------------------------------------------------------
library(shiny)
library(tidyverse)
library(ggsankey)
library(gridlayout)
library(bslib)
library(plotly)
library(DT)
library(ComplexHeatmap)
library(ggpubr)

source("helpers.R")
source("plotting.R")

# Reading data ------------------------------------------------------------


all_public_dss = readRDS("data/curated_all_public_metadata_nat_comm.rds")
training_dss = readRDS("data/Revision_updated_training_datasets.rds")
testing_dss = readRDS("data/Revision_updated_testing_datasets.rds")
all_eval_res = readRDS("data/all_eval_res.rds")
enrich_res = readRDS("data/context_enrich_results.rds")
HMDB_taxo_info = readRDS("data/hmdb_taxo_info.rds")
rel_scores = readRDS("data/Reliability_scores.rds")


# UI ----------------------------------------------------------------------
ui <- navbarPage(
  title = "METASPACE-ML: Context Explorer",
  selected = "Match context",
  collapsible = TRUE,
  theme = bslib::bs_theme(),
  tabPanel(
    title = "Match context",
    grid_container(
      layout = c(
        "ds_id match_context_1",
        "ds_id plotly_sunny   "
      ),
      row_sizes = c(
        "1fr",
        "1fr"
      ),
      col_sizes = c(
        "0.75fr",
        "1.25fr"
      ),
      gap_size = "10px",
      grid_card(
        area = "ds_id",
        card_body_fill(
          # radioButtons(
          #   inputId = "ds_id_check",
          #   label = "Do you have a METASPACE dataset ID ?",
          #   choices = list("Yes" = "y", "No" = "n"),
          #   width = "100%"
          # ),
          # conditionalPanel(
          #   condition = "input.ds_id_check == 'y'",
          #   textInput(
          #     inputId = "ds_id",
          #     label = "Dataset ID",
          #     value = "",
          #     width = "100%"
          #   )
          # ),
          card(
            full_screen = TRUE,
            card_header(span(strong("Select context parameters"))),
            card_body_fill(
              radioButtons(
                inputId = "KingdomId",
                label = "Kingdom",
                choices = list("Animal" = "Animal", "Plant" = "Plant"),
                width = "100%"
              ),
              radioButtons(
                inputId = "polarity",
                label = "Polarity",
                choices = list("Negative" = "Negative", "Positive" = "Positive"),
                width = "100%"
              ),
              selectInput(
                inputId = "SourceId",
                label = "Source",
                choices = list("MALDI" = "MALDI", "DESI" = "DESI",
                               "All" = "All")
              ),
              selectInput(
                inputId = "AnalyzerId",
                label = "Analyzer",
                choices = list("Orbitrap" = "Orbitrap", "FTICR" = "FTICR",
                               "All" = "All")
              ),
              conditionalPanel(
                condition = "input.KingdomId == 'Animal'",
                selectInput(
                  inputId = "animalspecies",
                  label = "Animal Species",
                  choices = list("Human" = "Homo sapiens",
                                 "Mouse" = "Mus musculus",
                                 "Other" = "OTHER"))
              ),
              conditionalPanel(
                condition = "input.KingdomId == 'Plant'",
                selectInput(
                  inputId = "planspecies",
                  label = "Plant Species",
                  choices = list(
                    "Populus" = "Populus",
                    "Sorghum" = "Sorghum",
                    "Other" = "OTHER"
                  ))
              ),
              numericInput(
                inputId = "minmzId",
                label = "min m/z",
                value = 5,
                min = 0,
                width = "100%"
              ),
              numericInput(
                inputId = "maxmzId",
                label = "max m/z",
                value = 1000,
                min = 50,
                width = "100%"
              ),
              selectInput(
                inputId = "sampletypeId",
                label = "Sample Type",
                choices = list("Tissue" = "Tissue", "Cells" = "Cells",
                               "Whole Organism" = "Whole Organism",
                               "Other" = "Uncurated",
                               "Uncertain" = "idk")
              )
            )
          ),
        )
      ),
      grid_card(
        area = "match_context_1",
        card_body_fill(
          textOutput(outputId = "context_desc"),
          plotOutput(outputId = "sankey", width = "100%")
        )
      ),
      grid_card(
        area = "plotly_sunny",
        card_body_fill(plotlyOutput(outputId = "sunny_plotly"))
      )
    )
  ),
  tabPanel(
    title = "Datasets",
    tabsetPanel(
      tabPanel(
        title = "Training",
        DTOutput(outputId = "trainingtableId", width = "95%"),
        downloadButton("downloadTrainingBtn", "Download")
      ),
      tabPanel(
        title = "Testing",
        DTOutput(outputId = "testingtableId", width = "95%"),
        downloadButton("downloadTestingBtn", "Download")
      )
    )
  ),
  tabPanel(
    title = "Evaluation",
    grid_container(
      layout = c(
        "area0 area1",
        "area2 area3"
      ),
      row_sizes = c(
        "1fr",
        "1fr"
      ),
      col_sizes = c(
        "1fr",
        "1fr"
      ),
      gap_size = "10px",
      grid_card(
        area = "area0",
        card_header(
          span(strong("Paired boxplot showing Area Under precision-recall Curve (AUC) per dataset. Check Figure 3 for more info. "))
        ),
        card_body_fill(plotOutput(outputId = "fig3a"))
      ),
      grid_card(
        area = "area1",
        card_header(
          strong("Boxplot showing the difference in number of annotations captured by METASPACE-ML relative  to MSM method per datasets")
        ),
        card_body_fill(plotOutput(outputId = "fig4a"))
      ),
      grid_card(
        area = "area2",
        card_header(
          strong("Relationship between the increased number of annotations compared to MSM and MAP scores for datasets in selected context")
        ),
        card_body_fill(plotOutput(outputId = "fig5a"))
      ),
      grid_card(
        area = "area3",
        card_header(strong("Reliability score distribution for optimal FDR thresholds. Check Figure 5D for more info.")),
        card_body_fill(
          plotOutput(outputId = "reliabilityId")
        )
      )
    )
  ),
  tabPanel(
    title = "Enrichment",
    grid_container(
      layout = c(
        "area0 area0",
        "area1 area2"
      ),
      row_sizes = c(
        "1fr",
        "1fr"
      ),
      col_sizes = c(
        "1.6fr",
        "0.4fr"
      ),
      gap_size = "10px",
      grid_card(
        area = "area0",
        card_body_fill(plotOutput(outputId = "fig7c", width = "100%",height = "100%"))
      ),
      grid_card(
        area = "area1",
        card_body_fill(
          DTOutput(outputId = "enrichtable", width = "95%")
        ),
        card_footer(
          downloadButton("downloadenrichBtn", "Download")
        )
      ),
      grid_card(
        area = "area2",
        card_header(strong("Metabolites in the selected dataset and enriched class")),
        card_body_fill(
          verbatimTextOutput("tp_markers")
          )
      )
    )
  )
)


# Server ------------------------------------------------------------------
server <- function(input, output, session) {

  kingdom = reactive(input$KingdomId)

  selected_organism <- reactive({
    switch(kingdom(),
           "Animal" = input$animalspecies,
           "Plant" = input$planspecies)
  })

  sel_contexts <- reactive({
    req(input$polarity, input$SourceId, input$AnalyzerId, input$minmzId,
        input$maxmzId, selected_organism(), input$sampletypeId)
    match_context(pol = input$polarity,
                  source = input$SourceId,
                  analyzer = input$AnalyzerId,
                  min_mz = input$minmzId,
                  max_mz = input$maxmzId,
                  organism = selected_organism(),
                  sample_type = input$sampletypeId
    )
  })

  training_context_df = reactive({
    selected_kingdom <- kingdom()
    training_ds_list <- training_dss[[selected_kingdom]][["30"]]
    get_context_df(ds_context_list = training_ds_list)
  })
  testing_context_df = reactive({
    selected_kingdom <- kingdom()
    testing_ds_list <- testing_dss[[selected_kingdom]]
    get_context_df(ds_context_list = testing_ds_list)
  })

  ds_id_exist = reactive(input$ds_id_check)

  sel_dss_training <- reactiveVal(NULL)
  sel_dss_testing <- reactiveVal(NULL)

  # observe({
  #   if(ds_id_exist() == "n"){
  #
  #   }
  #   else{
  #     #TODO write function to get sel context from api
  #   }
  #   print(length(sel_dss_testing()))
  # })

  observe({
    sel_contexts()

    filtered_training <- training_context_df()$ds_id[training_context_df()$context %in% sel_contexts()]
    filtered_testing <- testing_context_df()$ds_id[testing_context_df()$context %in% sel_contexts()]

    sel_dss_training(filtered_training)
    sel_dss_testing(filtered_testing)
  })

  train_table = reactive({

    req(sel_dss_training())

    all_public_dss %>%
      dplyr::filter(ds_id %in% sel_dss_training()) %>%
      dplyr::select(-is_public, -rp_range, -proj_id, -proj_name,
                    -submission_day, -regular_geom, -ibd_size) %>%
      add_hyperlink_ds()
  })

  testing_table = reactive({
    req(sel_dss_testing())

    all_public_dss %>%
      dplyr::filter(ds_id %in% sel_dss_testing()) %>%
      dplyr::select(-is_public, -rp_range, -proj_id, -proj_name,
                    -submission_day, -regular_geom, -ibd_size) %>%
      add_hyperlink_ds()
  })

  enrich_table = reactive({
    req(sel_contexts(), kingdom())
    prepare_enrich_table(enrich_res = enrich_res,
                           ds_list = testing_dss,
                           sel_context = sel_contexts(),
                           kingdom = kingdom(),
                           filter_by_adjpval = F,
                           min_TP = 3, use_FE = T)
  })

  output$trainingtableId = renderDT({

    train_table()

    # datatable(train_table,
    #           options = list(dom = 'Bfrtip', buttons = c('csv')),
    #           rownames = FALSE)

  }, escape = F)

  output$testingtableId = renderDT({

    testing_table()

    # datatable(testing_table,
    #           options = list(dom = 'Bfrtip', buttons = c('csv')),
    #           rownames = FALSE)

  }, escape = F)

  output$downloadTrainingBtn <- downloadHandler(
    filename = function() {
      paste("training_data", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(train_table(), file)
    }
  )

  output$downloadTestingBtn <- downloadHandler(
    filename = function() {
      paste("testing_data", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(testing_table(), file)
    }
  )

  output$sankey <- renderPlot({

    req(testing_context_df(),kingdom(),sel_dss_testing())

    plot_new_sankey(meta_df = testing_context_df(),
                    kingdom = kingdom(),
                    sel_dss = sel_dss_testing(),
                    plot_title = "Testing datasets")
  })

  output$sunny_plotly <- renderPlotly({
    req(testing_context_df(),sel_dss_testing())

    plot_sunburst_meta(meta_df = testing_context_df(),
                       sel_dss = sel_dss_testing())

  })

  output$fig3a <- renderPlot({

    req(kingdom(), sel_contexts(),testing_context_df())

    if(length(intersect(testing_context_df()$context, sel_contexts())) == 1){
      Performance_plots_pipeline(eval_res = all_eval_res,
                                 ds_list = testing_dss,
                                 sel_context = sel_contexts(),
                                 data_type = "Testing",
                                 context_size = 30,
                                 kingdom = kingdom(),
                                 mean_eval_context = F,
                                 metrics_per_grp = F,
                                 FDR.pct = NULL,
                                 plot_type = "paired_box",
                                 x_column = "score_type",
                                 y_column = "metric_value",
                                 color_column = "score_type",
                                 hide_x_axis = T,hide_x_axis_label = T,
                                 context_specific = F,
                                 in_plot_text_size = 3,
                                 plot_title = paste0("Testing datasets", "\n", "MAP"))
    }
    else if (length(intersect(testing_context_df()$context, sel_contexts())) > 1){

      AP_box_animal = reactive({
        prepare_plot_df(eval_res = all_eval_res,
                        ds_list = testing_dss,
                        sel_context = sel_contexts(),
                        data_type = "Testing",context_size = 30,
                        kingdom = kingdom(),metrics_per_grp = F,
                        mean_eval_context = F)
      })

      make_density_hm_context(df = AP_box_animal(), FDR = 10,
                              vals = "MAP", kingdom = kingdom(),
                              plot_title = paste0("Testing datasets", "\n", "MAP"))
    }
  })

  output$fig4a <- renderPlot({

    req(kingdom(), sel_contexts(),testing_context_df())

    if(length(intersect(testing_context_df()$context, sel_contexts())) == 1){
      Performance_plots_pipeline(eval_res = all_eval_res,
                                 ds_list = testing_dss,
                                 sel_context = sel_contexts(),
                                 data_type = "Testing",
                                 context_size = 30,
                                 kingdom = kingdom(),
                                 mean_eval_context = F,
                                 metrics_per_grp = F,
                                 FDR.pct = NULL,
                                 plot_type = "box",
                                 x_column = "FDR_pct",
                                 y_column = "LogDiff",
                                 color_column = "FDR_pct",
                                 hide_x_axis = F,hide_x_axis_label = F,
                                 context_specific = F,
                                 in_plot_text_size = 5,
                                 plot_title = paste0("Testing datasets", "\n", "Annotations"))
    }
    else if (length(intersect(testing_context_df()$context, sel_contexts())) > 1){
      Annot_MAP_animal = reactive({
        prepare_plot_df(eval_res = all_eval_res,
                        ds_list = testing_dss,
                        sel_context = sel_contexts(),
                        data_type = "Testing",context_size = 30,
                        kingdom = kingdom(),metrics_per_grp = F,
                        mean_eval_context = F,FDR.pct = 10)
      })

      make_density_hm_context(df = Annot_MAP_animal(),FDR = 10,vals = "LogDiff",
                              kingdom = kingdom(), plot_title = NULL)
    }
  })

  output$fig5a <- renderPlot({

    req(kingdom(), sel_contexts(),testing_context_df())

    Performance_plots_pipeline(eval_res = all_eval_res,
                               ds_list = testing_dss,
                               sel_context = sel_contexts(),
                               data_type = "Testing",
                               context_size = 30,
                               kingdom = kingdom(),
                               mean_eval_context = F,
                               metrics_per_grp = F,
                               FDR.pct = 10,
                               plot_type = "2d_bin",
                               x_column = "LogDiff",
                               y_column = "metric_value",
                               color_column = "score_type",
                               hide_x_axis = F,hide_x_axis_label = F,
                               context_specific = F,
                               in_plot_text_size = 5)
  })

  output$reliabilityId <- renderPlot({

    req(kingdom(), sel_dss_testing())

    ggstats_wrapper(df = rel_scores,
                    kingdom = kingdom(),
                    sel_dss = sel_dss_testing(),
                    x_column = "FDR_pct",
                    y_column = "rel_score",paired = F,
                    centrality.plotting = F,
                    results.subtitle = F, parametric = F) +
      ylab("Reliability Score") +
      xlab("FDR (%)") +
      theme(plot.title = element_text(hjust = 0.5),
            axis.title = element_text(size = 14)) +
      ylim(c(0,1))
  })

  output$fig7c <- renderPlot({
    req(kingdom(), sel_contexts(),testing_context_df(),sel_dss_testing())

    if(length(intersect(testing_context_df()$context, sel_contexts())) == 1){

      plot_enrich_boxplot(enrich_res = enrich_res,
                          kingdom = kingdom(),
                          sel_dss = sel_dss_testing(),
                          HMDB_taxo_info = HMDB_taxo_info,
                          pval_thresh = 0.05, min_TP = 3,
                          min_dss_prop = 0.1)

    }
    else if (length(intersect(testing_context_df()$context, sel_contexts())) > 1){

      plot_enrich_hm_context(enrich_res = enrich_res,
                             ds_list = testing_dss,
                             sel_context = sel_contexts(),
                             kingdom = kingdom(),
                             filter_by_adjpval = F,
                             min_TP = 3, use_FE = T,
                             min_dss_prop = 0.1)
    }
  })

  output$enrichtable <- renderDT({
    enrich_table() %>%
      dplyr::select(-TP_markers)
  },selection = "single")

  output$downloadenrichBtn <- downloadHandler(
    filename = function() {
      paste("enrichment_data", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(enrich_table(), file)
    }
  )

  output$tp_markers <- renderPrint({
    enrich_table()[input$enrichtable_rows_selected, "TP_markers"] %>%
      as.character() %>%
      strsplit(",") %>%
      unlist() %>%
      paste(sep = "", collapse = "\n") %>%
      cat()
  })


}

shinyApp(ui, server)


