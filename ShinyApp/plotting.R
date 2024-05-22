
# General -----------------------------------------------------------------

get_d_hm_colors_context = function(kingdom = c("Animal", "Plant"), add_mz_class = F){
  if (kingdom == "Animal"){
    cols = list("polarity" = c("Positive" = "#A6CEE3", "Negative" = "#1F78B4"),
                "analyzer" = c("Orbitrap" = "#B2DF8A", "FTICR" = "#33A02C"),
                "source" = c("MALDI" = "#FB9A99", "DESI" = "#E31A1C"),
                "Species" = c("Homo sapiens" = "#8DD3C7", "Mus musculus" = "#FFFFB3", "OTHER" = "#BC80BD"),
                "Sample_type" = c("Tissue" = "#FDBF6F", "Uncurated" = "#CAB2D6", "Whole Organism" = "#DFC497", "Cells" = "#664098"))
  }
  else{
    cols = list("polarity" = c("Positive" = "#A6CEE3", "Negative" = "#1F78B4"),
                "analyzer" = c("Orbitrap" = "#B2DF8A", "FTICR" = "#33A02C"),
                "source" = c("MALDI" = "#FB9A99", "DESI" = "#E31A1C"),
                "Species" = c("Populus" = "#8DD3C7", "Sorghum" = "#FFFFB3", "OTHER" = "#BC80BD"),
                "Sample_type" = c("Tissue" = "#FDBF6F", "Uncurated" = "#CAB2D6"))
  }
  if (add_mz_class){
    cols[["mz_class"]] = c("Lipids" = "#D1EAC6", "Lipids_and_small_moleclues" = "#62B5C2")
  }
  
  return(cols)
}
sigFunc <- function(x) {
  ifelse(x < 0.001, "***",
         ifelse(x < 0.01, "**",
                ifelse(x < 0.05, "*", "ns")))
}


# Descriptive -------------------------------------------------------------


plot_new_sankey <- function(meta_df,ds_context_list = NULL,
                            kingdom, sel_dss, plot_title) {
  
  if(!is.null(ds_context_list)) {
    context_df = ds_context_list %>%
      unlist() %>% as.data.frame() %>%
      rownames_to_column("context")
    
    colnames(context_df)[2] = "ds_id"
    context_df$context = sub(".*[.]", "", context_df$context)
    context_df$context = sub("[0-9]+.*", "", context_df$context)
    
    context_df = context_df %>%
      distinct() %>%
      tidyr::separate(col = context,
                      into = c("polarity","mz_class",
                               "analyzer", "source",
                               "Species", "Sample_type"),
                      sep = ",", remove = T)
  }
  else{
    context_df = meta_df
  }
  
  fig <- sankey_plot_dss(
    meta_df = context_df,
    sel_dss = sel_dss,
    plot_title = plot_title,
    kingdom = kingdom
  )
  return(fig)
}

sankey_plot_dss = function(meta_df, sel_dss, plot_title, kingdom){
  
  if ("mz_class" %in% colnames(meta_df)) {
    if (length(unique(meta_df$mz_class)) > 1){
      cols = get_d_hm_colors_context(kingdom = kingdom,add_mz_class = T) %>% unlist()
      context_cols = c("polarity","mz_class","source","analyzer","Species","Sample_type")
      meta_df = meta_df[which(meta_df$ds_id %in% sel_dss),context_cols]
      meta_df_long = meta_df %>% ggsankey::make_long(polarity,mz_class, source,
                                                     analyzer,Species,Sample_type)
    }
    else{
      cols = get_d_hm_colors_context(kingdom = kingdom,add_mz_class = F) %>% unlist()
      context_cols = c("polarity","source","analyzer","Species","Sample_type")
      meta_df = meta_df[which(meta_df$ds_id %in% sel_dss),context_cols]
      meta_df_long = meta_df %>% ggsankey::make_long(polarity, source,
                                                     analyzer,Species,Sample_type)
    }
  }
  else{
    cols = get_d_hm_colors_context(kingdom = kingdom,add_mz_class = F) %>% unlist()
    context_cols = c("polarity","source","analyzer","Species","Sample_type")
    meta_df = meta_df[which(meta_df$ds_id %in% sel_dss),context_cols]
    meta_df_long = meta_df %>% ggsankey::make_long(polarity, source,
                                                   analyzer,Species,Sample_type)
  }
  names(cols) = sub(".*[.]","",names(cols))
  
  p = ggplot(meta_df_long, aes(x = x, next_x = next_x, node = node, next_node = next_node,
                               fill = factor(node), label = node)) +
    geom_alluvial(flow.alpha = .6, space = 1) +
    geom_alluvial_text(color = "Black", angle = 0) +
    scale_fill_manual(values = cols) +
    #scale_fill_brewer(palette = "Paired", type = "qual") +
    theme_alluvial(base_size = 15) +
    labs(x = NULL) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = .5),
          axis.text.x = element_text(size = 15, color = "Black"),
          axis.title.y = element_text(size = 15, color = "Black"),
          axis.text.y = element_text(size = 15, color = "Black")) +
    # scale_x_discrete(labels = addline_format(colnames(train_ds))) +
    ggtitle(plot_title) +
    ylab("Number of Datasets")
  
  return(p)
  
}

plot_sunburst_meta <- function(meta_df, sel_dss = NULL) {
  context_cols = c("polarity","source","analyzer","Species","Sample_type")
  if(!is.null(sel_dss)){
    meta_df = meta_df[meta_df$ds_id %in% sel_dss,context_cols]
  }
  
  trial <- meta_df %>%
    dplyr::count(polarity, source, analyzer, Sample_type, Species) %>%
    dplyr::distinct() %>%
    plotme:::create_all_col_params(F, T)
  
  trial$labels <- gsub("_", "<br>", trial$labels)
  trial$labels <- paste0(trial$labels, "\n", trial$values)
  
  
  
  x <- purrr::exec(plotly::plot_ly, !!!trial,
                   type = "sunburst",
                   branchvalues = "total"
  )
  return(x)
}


# Evaluation --------------------------------------------------------------
get_ggmat_genbase <- function(df, plot_type = c(
  "box", "bar", "scatter",
  "intens_box", "2d_bin"
),
x_column = "score_type",
y_column = "metric_value",
color_column = "score_type", in_plot_text_size = 2) {
  
  if(x_column == "FDR_pct" & plot_type == "box") {
    df$FDR_pct = factor(df$FDR_pct, levels = c("5", "10", "20", "50"))
    if (y_column %in% c("LogDiff", "LFC")){
      df = df[df$score_type == "METASPACE-ML",]
    }
  }
  if(plot_type == "scatter" | plot_type == "2d_bin"){
    df = df[df$score_type == "METASPACE-ML",]
  }
  
  if (plot_type == "box") {
    gg_base <- ggplot(
      data = df,
      aes(
        x = .data[[x_column]], y = .data[[y_column]],
        color = .data[[color_column]]
      )
    ) +
      geom_boxplot(
        width = 0.5, position = "dodge",
        notch = F, outlier.color = "Black", outlier.shape = NA) +
      geom_jitter(size = 0.5, alpha = 0.3, width = 0.1) +
      stat_summary(fun=mean, geom="errorbar", aes(ymax = ..y..,
                                                  ymin = ..y..), size = 1,
                   linetype = "dashed", width = 0.5)
  } else if (plot_type == "bar") {
    gg_base <- ggplot(
      data = df,
      aes(
        x = .data[[x_column]], y = .data[[y_column]],
        fill = .data[[color_column]]
      )
    ) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_text(aes(label = round(.data[[y_column]], 2)),
                position = position_dodge(width = 0.9),
                vjust = 1, hjust = 0, size = in_plot_text_size
      ) +
      geom_errorbar(aes(ymin=.data[[y_column]]-.data[["metric_sd"]], ymax=.data[[y_column]]+.data[["metric_sd"]]),
                    width=.2,
                    position=position_dodge(.9)) +
      scale_fill_manual(values = c("#A6CEE3", "#FB9A99")) +
      theme(legend.position = "none")
  } else if (plot_type == "scatter") {
    gg_base <- ggplot(
      data = df,
      aes(
        x = .data[[x_column]], y = .data[[y_column]],
        color = .data[[color_column]]
      )
    ) +
      geom_point(size = 2, alpha = 0.6) +
      # geom_line(aes(group=.data[[color_column]]), linetype=2) +
      # geom_smooth(method = "lm", se = F, size = 1) +
      # scale_color_manual(values = c("#A6CEE3", "#FB9A99")) +
      theme(legend.position = "none")
  } else if (plot_type == "intens_box") {
    gg_base <- ggplot(data = df, mapping = aes(
      x = .data[[x_column]], y = .data[[y_column]],
      color = log10(.data[[color_column]])
    )) +
      geom_boxplot(outlier.shape = NA, size = 1) +
      geom_jitter(
        size = 3, alpha = 0.8,
        position = position_jitterdodge(dodge.width = 0.9)
      ) +
      scale_y_continuous(n.breaks = 10, limits = c(0, 10)) +
      scale_color_continuous(type = "viridis", name = "Log10 (Number of ions)") +
      # scale_size_binned_area(breaks = c(1,10,50,100,300,500),
      #                        name = "Number of ions", max_size = 10) +
      theme_pubr() +
      # theme(legend.text = element_text(size = 14),
      #       legend.title = element_text(size = 16),
      #       axis.title = element_text(size = 16),
      #       axis.text = element_text(size = 14))+
      geom_signif(
        comparisons = list(c("MSM", "METASPACE-ML_only")),
        map_signif_level = sigFunc, na.rm = T, textsize = 8,
        y_position = c(10, 8, 9, 10)
      )
  } else if (plot_type == "2d_bin"){
    gg_base = ggplot(df, aes(x = .data[[x_column]], y = .data[[y_column]])) +
      geom_bin_2d(binwidth = c(0.5, 0.1)) +
      stat_bin2d(geom = "text", aes(label = ..count..),
                 binwidth = c(0.5, 0.1), size = in_plot_text_size) +
      scale_fill_gradient(low = "snow", high = "red") +
      xlim(-3, 3) +
      ylim(0, 1) +
      scale_y_continuous(breaks = seq(0, 1, 0.2)) +
      scale_x_continuous(breaks = seq(-3, 3, 1)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none") +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted")
  }
  else {
    stop("Invalid plot type")
  }
  return(gg_base)
}

Performance_plots_pipeline = function(eval_res,
                                      ds_list,
                                      sel_context,
                                      data_type = c("Training", "Testing"),
                                      context_size = 30,
                                      kingdom = c("Animal", "Plant"),
                                      mean_eval_context = F,
                                      metrics_per_grp = F,
                                      FDR.pct = NULL,
                                      plot_type = c("box", "bar", "scatter",
                                                    "intens_box", "paired_box"),
                                      x_column = "score_type",
                                      y_column = "metric_value",
                                      color_column = "score_type",
                                      hide_x_axis = F,hide_x_axis_label = T,
                                      context_specific = T,
                                      in_plot_text_size = 2,
                                      plot_title = NULL){
  
  eval_df = eval_res[[kingdom]][["30"]][[data_type]][["eval_df"]]
  annot_df = eval_res[[kingdom]][["30"]][[data_type]][["annot_df"]]
  
  y_label = switch(y_column, "metric_value" = "MAP", "FDR_pct" = "FDR (%)",
                   "Difference" = "Difference", "LogDiff" = "(+-)Log10(Absolute Difference)",
                   "LFC" = "Log2 Fold Change ")
  x_label = switch(x_column, "metric_value" = "MAP", "FDR_pct" = "FDR (%)",
                   "Difference" = "Difference", "LogDiff" = "(+-)Log10(Absolute Difference)",
                   "LFC" = "Log2 Fold Change", "score_type" = "Score type")
  
  if(!context_specific){
    plot_df = prepare_plot_df(eval_res = eval_res,ds_list = ds_list,
                              sel_context = sel_context,
                              data_type = data_type, context_size = context_size,
                              kingdom = match.arg(kingdom),metrics_per_grp = metrics_per_grp,
                              mean_eval_context = F,FDR.pct = FDR.pct)
    plot_df = plot_df[,!colnames(plot_df) %in% c("context","polarity","mz_class",
                                                 "analyzer", "source",
                                                 "Species", "Sample_type", "source_analyzer")] %>%
      dplyr::distinct() %>%
      dplyr::filter(!is.na(.data[[x_column]]), !is.na(.data[[y_column]]))
    if (plot_type == "bar"){
      plot_df = plot_df %>%
        dplyr::group_by(.data[[x_column]], metric_type) %>%
        dplyr::summarise(metric_sd = sd(metric_value,na.rm = T),
                         metric_value = mean(metric_value, na.rm = T)) %>%
        dplyr::ungroup()
    }
    
    if (plot_type == "paired_box"){
      PR_Data = plot_df %>%
        dplyr::select(ds_id,{{x_column}},{{y_column}}) %>%
        dplyr::distinct() %>%
        tidyr::spread(key = .data[[x_column]], value = .data[[y_column]])
      
      PR_Data$Delta = ifelse(PR_Data$`METASPACE-ML` > PR_Data$MSM, "Higher", "Lower")
      PR_Data = PR_Data[!is.na(PR_Data$Delta),]
      
      Delta_count = PR_Data %>%
        dplyr::group_by(Delta) %>%
        dplyr::summarise(n_ds = n_distinct(ds_id))
      
      PR_Data = PR_Data %>%
        dplyr::left_join(Delta_count, by = "Delta") %>%
        dplyr::mutate(Delta = paste(Delta, "\n(", n_ds, ")", sep = ""))
      
      p = PR_Data %>% ggpubr::ggpaired(cond1 = "METASPACE-ML", cond2 = "MSM",
                               color = "condition", linetype = "dotted",
                               line.size = 0.5, point.size = 2,
                               palette = c("#A6CEE3", "#FB9A99"),
                               ggtheme = theme_pubr(legend = "bottom"),
                               facet.by = "Delta",line.color = "#F4C08B",
                               xlab = FALSE, ylab = "PR-AUC", title = "Paired PR-AUC per dataset",
                               short.panel.labs = T) +
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(size = 20),
              axis.text.y = element_text(size = 20),
              axis.title.y = element_text(size = 20),
              legend.text = element_text(size = 20),
              legend.title = element_blank(),
              strip.text = element_text(size = 20),
              strip.background = element_rect(fill = "white")) +
        ggpubr::stat_compare_means(aes(group = condition), label = "p.signif" ,paired = F,
                           size = 8, label.x.npc = "center", label.y.npc = 0.95)
    }
    else{
      gg_base = get_ggmat_genbase(df = plot_df, plot_type = plot_type, x_column = x_column,
                                  y_column = y_column, color_column = color_column,
                                  in_plot_text_size = in_plot_text_size)
      
      p = gg_base +
        ggpubr::theme_pubr() +
        ylab(y_label) +
        xlab(x_label) +
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 14),
              axis.title.x = element_text(size = 14))
    }
  }
  else{
    plot_df = prepare_plot_df(eval_df = eval_df,
                              annot_df = annot_df,
                              data_type = data_type, context_size = 30,
                              kingdom = match.arg(kingdom),metrics_per_grp = metrics_per_grp,
                              mean_eval_context = mean_eval_context,FDR.pct = FDR.pct) %>%
      filter(!is.na(.data[[x_column]]), !is.na(.data[[y_column]]))
    
    p = ggmat_wrapper(plot_df = plot_df, plot_type = plot_type,
                      x_column = x_column, y_column = y_column,
                      color_column = color_column ,hide_x_axis = hide_x_axis,
                      hide_x_axis_label = hide_x_axis_label,
                      in_plot_text_size = in_plot_text_size)
  }
  if(hide_x_axis){
    p = p + theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank())
  }
  if(hide_x_axis_label){
    p = p + theme(axis.title.x = element_blank())
  }
  if(!is.null(plot_title)){
    p = p + ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5))
  }
  return(p)
}


make_density_hm_context = function(df, FDR = 10,
                                   vals = c("Diff", "LogDiff", "LFC", "MAP"),
                                   kingdom = c("Animal", "Plant"), col_split = NULL,
                                   plot_title){
  
  val_col = switch(vals, "Diff" = "Difference", "LogDiff" = "LogDiff",
                   "LFC" = "LFC", "MAP" = "metric_value")
  
  anno_df = df %>%
    dplyr::select(context) %>%
    as.data.frame(check.names = F) %>%
    tidyr::separate(col = context,
                    into = c("polarity","mz_class",
                             "analyzer", "source",
                             "Species", "Sample_type"),
                    sep = ",", remove = F) %>%
    dplyr::select(-mz_class) %>%
    dplyr::distinct() %>%
    tibble::column_to_rownames(var = "context")
  
  ds_context = df %>%
    dplyr::select(ds_id, context) %>%
    dplyr::distinct()
  
  d_hm_mat = df %>%
    dplyr::filter(FDR_pct == FDR,
                  ds_id %in% ds_context$ds_id,
                  score_type == "METASPACE-ML")
  d_hm_mat = d_hm_mat[,which(colnames(d_hm_mat) %in% c("ds_id", "context", val_col))]
  colnames(d_hm_mat)[which(colnames(d_hm_mat) == val_col)] = "values"
  d_hm_mat = d_hm_mat %>%
    tidyr::pivot_wider(names_from = context, values_from = values) %>%
    tibble::column_to_rownames(var = "ds_id") %>%
    as.matrix()
  
  dm_cols = get_d_hm_colors_context(kingdom = kingdom)
  
  
  p =
    ComplexHeatmap::densityHeatmap(data = d_hm_mat, show_column_names = F, cluster_columns = T,
                   show_quantiles = T, title = plot_title,
                   ylab = vals, quantile_gp = gpar(fontsize = 10),
                   show_heatmap_legend = T,column_split = col_split
    ) %v%
    ComplexHeatmap::HeatmapAnnotation(df = anno_df,
                      col = dm_cols,show_legend = T,
                      annotation_legend_param = list("labels_gp" =
                                                       gpar(fontsize = 12),
                                                     "title_gp" =
                                                       gpar(fontsize = 14, fontface = "bold")))
  
  
  return(p)
  
  
}

ggstats_wrapper = function(df,
                           kingdom,
                           sel_dss,
                           x_column, y_column,
                           paired = F, centrality.plotting = T,
                           results.subtitle = F,parametric = F,
                           sig_only = F){
  
  df = df[[kingdom]] %>% 
    dplyr::filter(ds_id %in% sel_dss)
  
  y_label = switch(y_column, "metric_value" = "MAP", "FDR_pct" = "FDR (%)",
                   "Difference" = "Difference", "LogDiff" = "(+-)Log10(Absolute Difference)",
                   "LFC" = "Log2 Fold Change ")
  x_label = switch(x_column, "metric_value" = "MAP", "FDR_pct" = "FDR (%)",
                   "Difference" = "Difference", "LogDiff" = "(+-)Log10(Absolute Difference)",
                   "LFC" = "Log2 Fold Change", "score_type" = "Score type")
  
  ggstats_func = ifelse(paired, ggstatsplot::ggwithinstats, ggstatsplot::ggbetweenstats)
  
  param_type = ifelse(parametric, "parametric", "nonparametric")
  
  p = ggstats_func(data = df,x = {{x_column}},y = {{y_column}},type = param_type,
                   p.adjust.method = "BH", pairwise.display = "none",
                   results.subtitle = results.subtitle ,pairwise.comparisons = T,
                   centrality.plotting = centrality.plotting)
  if(length(unique(df[[x_column]])) == 2){
    return(p)
  }
  
  p_comparison = ggstatsplot::extract_stats(p)[["pairwise_comparisons_data"]] %>%
    dplyr::mutate(groups = purrr::map2(group1, group2, function(x,y){c(x,y)})) %>%
    dplyr::mutate(p_astersik = sigFunc(p.value))
  
  if (sig_only){
    p_comparison = p_comparison %>%
      dplyr::filter(p_astersik != "ns")
  }
  
  g = ggstats_func(data = df,x = {{x_column}},y = {{y_column}},type = param_type,
                   p.adjust.method = "BH", pairwise.display = "none",
                   results.subtitle = F,pairwise.comparisons = F,
                   centrality.plotting = T, plot.type = "box")
  
  final_p = g + ggsignif::geom_signif(
    comparisons = p_comparison$groups,
    map_signif_level = T,
    annotations = p_comparison$p_astersik,
    test = NULL,
    y_position = ggstatsplot:::.ggsignif_xy(pull(df,.data[[x_column]])
                                            , pull(df,.data[[y_column]])),
    na.rm = T,
    textsize = 5
  ) +
    scale_y_continuous(n.breaks = 10) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 12)) +
    xlab(x_label) +
    ylab(y_label)
  
  return(final_p)
  
}


# Enrichment --------------------------------------------------------------
plot_enrich_boxplot = function(enrich_res,sel_dss,
                               kingdom = c("Animal", "Plant"),
                               HMDB_taxo_info,
                               pval_thresh = 0.05, min_TP = 3,
                               min_dss_prop = 0.1){
  plot_data = enrich_res[[kingdom]] %>% 
    dplyr::filter(pval < pval_thresh, TP >= min_TP, ds_id %in% sel_dss)
  n_sig_dss = length(unique(plot_data$ds_id))
  
  plot_data$OR = log2(plot_data$FE)
  
  
  filtered = plot_data %>% 
    dplyr::group_by(term) %>%
    dplyr::summarise(n_ds = n(), avg_LOR = mean(OR))
  
  filtered = filtered[which(filtered$n_ds >= ceiling(min_dss_prop * n_sig_dss)),]
  filtered = filtered[which(filtered$term != ""),]
  
  
  x = HMDB_taxo_info[which(HMDB_taxo_info$sub_class %in% filtered$term),]
  x = x %>% 
    dplyr::select(class, sub_class) %>% 
    dplyr::distinct()
  filtered = filtered %>% left_join(x, by = c("term" = "sub_class"))
  
  plot_data = plot_data[which(plot_data$term %in% filtered$term),]
  plot_data = plot_data %>%
    dplyr::left_join(filtered) %>%
    dplyr::select(term, class, ds_id, n_ds, OR)
  
  plot_data = plot_data[order(plot_data$class),]
  
  cols = RColorBrewer::brewer.pal(length(unique(plot_data$class)),"Paired")
  
  plot_data$myaxis = paste0(plot_data$term, " (", plot_data$n_ds, "/",
                            n_sig_dss,")")
  group_ordered = with(plot_data, reorder(myaxis, OR, median))
  
  plot_data$myaxis = factor(plot_data$myaxis, levels = levels(group_ordered))
  
  p = plot_data %>%
    ggplot(aes(x=myaxis, y=OR, fill=class)) +
    # geom_violin(width=1.4) +
    geom_boxplot(width= 1, alpha=0.6) +
    # geom_jitter(size= 0.5, alpha=0.5) +
    scale_fill_manual(values= cols) +
    scale_y_continuous(n.breaks = 16) +
    ggpubr::theme_pubr(legend = "bottom") +
    theme(plot.title = element_text(size=18, hjust = 0.5),
          plot.subtitle = element_text(size=16, hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size=18, hjust = 0.5),
          legend.text = element_text(size=18),
          legend.title = element_text(size = 16),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18)) +
    coord_flip() +
    geom_hline(yintercept =  0, linetype="dashed", color = "red") +
    ggtitle("Enrichment of Metabolite classes",
            subtitle = "Annotations only captured by METASPACE-ML") +
    ylab("Log2 Fold Enrichment")
  
  return(p)
}

plot_enrich_hm_context = function(enrich_res, ds_list,
                                  sel_context,
                                  filter_by_adjpval = T,
                                  min_TP = 3, use_FE = T,
                                  kingdom = c("Animal", "Plant"),
                                  min_dss_prop = 0.1){
  enrich_res = enrich_res[[kingdom]] %>%
    dplyr::filter(TP > min_TP, term != "") %>%
    dplyr::group_by(ds_id) %>%
    dplyr::mutate(adj_p_val = p.adjust(pval, method = "BH"))
  
  if (filter_by_adjpval){
    enrich_res = enrich_res %>%
      dplyr::filter(adj_p_val < 0.05)
  }
  else{
    enrich_res = enrich_res %>%
      dplyr::filter(pval < 0.05)
  }
  if (use_FE){
    enrich_res$OR = log2(enrich_res$FE)
  }
  else{
    enrich_res$OR = log2(enrich_res$OR)
  }
  
  ds_context_list = ds_list[[match.arg(kingdom)]]
  names(ds_context_list) = stringr::str_replace_all(names(ds_context_list),
                                                    c("ESI" = "DESI"))
  
  context_df = ds_context_list %>%
    as.data.frame(check.names = F) %>%
    gather(key = "context", value = "ds_id") %>%
    tidyr::separate(col = context,
                    into = c("polarity","mz_class",
                             "analyzer", "source",
                             "Species", "Sample_type"),
                    sep = ",", remove = F) %>% 
    dplyr::filter(context %in% sel_context)
  
  
  plot_df = enrich_res %>% dplyr::left_join(context_df, by = "ds_id") %>%
    dplyr::select(ds_id, context, term, OR) %>%
    dplyr::distinct() %>%
    dplyr::group_by(term, context) %>%
    dplyr::summarise(OR = mean(OR),
                     n_sig_ds = n_distinct(ds_id),.groups = "drop") %>%
    dplyr::filter(!is.na(context))
  
  wide_mat = plot_df %>%
    dplyr::filter(n_sig_ds >= min_dss_prop * 30) %>%
    dplyr::select(term, context, OR) %>%
    pivot_wider(names_from = context, values_from = OR) %>%
    column_to_rownames(var = "term") %>%
    as.matrix()
  
  anno_df = plot_df %>%
    dplyr::filter(n_sig_ds >= min_dss_prop * 30) %>% 
    dplyr::select(context) %>%
    dplyr::distinct() %>%
    tidyr::separate(col = context,
                    into = c("polarity","mz_class",
                             "analyzer", "source",
                             "Species", "Sample_type"),
                    sep = ",", remove = F) %>%
    dplyr::filter(!is.na(context)) %>% 
    tibble::column_to_rownames(var = "context")
  
  label_mat = plot_df %>%
    dplyr::filter(n_sig_ds >= min_dss_prop * 30) %>%
    dplyr::select(term, context, n_sig_ds) %>%
    tidyr::pivot_wider(names_from = context, values_from = n_sig_ds) %>%
    tibble::column_to_rownames(var = "term") %>%
    as.matrix()
  
  label_mat[is.na(label_mat)] = ""
  
  colors = get_d_hm_colors_context(kingdom = match.arg(kingdom),
                                   add_mz_class = T)
  
  wide_mat[is.na(wide_mat)] = 0
  
  ComplexHeatmap::pheatmap(mat = wide_mat,
                           color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu")))(100),
                           display_numbers = label_mat,annotation_col = anno_df,
                           show_colnames = F,cluster_rows = F, cluster_cols = F,
                           annotation_colors = colors,number_color = "black",
                           legend = T,annotation_legend = T,
                           fontsize_number = 12,fontsize = 12,
                           main = paste0("Enrichment per context\n", match.arg(kingdom)),
                           na_col = "white", border_color = "white",
                           heatmap_legend_param = list(title = "Log2 fold enrichment"))
  
}