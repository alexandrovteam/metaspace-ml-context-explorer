get_mz_class = function(min_mz, max_mz){
  if(max_mz <= 400){
    return("Small_molecules")
  }
  if(min_mz > 500){
    return("Lipids")
  }
  return("Lipids_and_small_moleclues")
}

match_context = function(pol, source, analyzer, 
                         min_mz, max_mz, organism,
                         sample_type){
  
  source = stringr::str_replace(source, "DESI", "ESI")
  
  if(source == "All"){source = c("MALDI", "ESI")}
  source = stringr::str_replace_all(source,
                                    c("ESI" = "DESI"))
  if(analyzer == "All"){analyzer = c("Orbitrap", "FTICR")}
  mz_class = get_mz_class(min_mz = min_mz, max_mz = max_mz)
  
  if(sample_type == "idk"){sample_type = c("Tissue", "Cells","Whole Organism",
                                           "Uncurated")}
  
  sel_contexts = list(
    polarity = pol,
    class = mz_class, 
    analyz = analyzer,
    ion_soruce = source,
    org = organism,
    type = sample_type
  ) %>% 
    expand.grid() %>% 
    apply(1, paste, collapse = ",")
  
  return(sel_contexts)
  
}

get_context_df = function(ds_context_list){
  
  names(ds_context_list) = stringr::str_replace_all(names(ds_context_list),
                                                   c("ESI" = "DESI"))
  context_df = ds_context_list %>%
    as.data.frame(check.names = F) %>%
    tidyr::gather(key = "context", value = "ds_id") %>%
    tidyr::separate(col = context,
                    into = c("polarity","mz_class",
                             "analyzer", "source",
                             "Species", "Sample_type"),
                    sep = ",", remove = F)
  return(context_df)
}

add_hyperlink_ds = function(df){
  df$url = paste0("https://metaspace2020.eu/dataset/", df$ds_id)
  df$url <- paste0("<a href='",df$url,"'>",df$url,"</a>")
  df %>% 
    dplyr::relocate(url, .after = ds_id)
}

prepare_plot_df = function(eval_res, ds_list,
                           sel_context,
                           data_type = c("Training", "Testing"),
                           context_size = 30,
                           kingdom = c("Animal", "Plant"),
                           mean_eval_context = F,
                           metrics_per_grp = F,
                           FDR.pct = NULL){
  eval_df = eval_res[[kingdom]][["30"]][[data_type]][["eval_df"]]
  annot_df = eval_res[[kingdom]][["30"]][[data_type]][["annot_df"]]
  
  if (data_type == "Training"){
    ds_context_list = ds_list[[match.arg(kingdom)]][[as.character(context_size)]]
  } else if (data_type == "Testing"){
    ds_context_list = ds_list[[match.arg(kingdom)]]
  }
  
  context_df = get_context_df(ds_context_list)
  context_df = context_df[context_df$context %in% sel_context,]
  
  eval_df$score_type = str_replace_all(eval_df$score_type,
                                       c("Catboost" = "METASPACE-ML"))
  eval_df = eval_df %>%
    tidyr::separate(group_name, into = c("ds_id", "adduct"), sep = ",",
                    remove = F) %>%
    dplyr::filter(ds_id %in% context_df$ds_id) %>% 
    dplyr::left_join(context_df, by = "ds_id")
  
  if(!metrics_per_grp){
    eval_df = eval_df %>%
      dplyr::group_by(context, score_type, metric_type, ds_id) %>%
      dplyr::summarise(metric_sd = sd(metric_value,
                                      na.rm = T),
                       metric_value = mean(metric_value,na.rm = T)) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(context_df, by = c("ds_id", "context"))
  }
  
  if(mean_eval_context){
    eval_df = eval_df %>%
      dplyr::group_by(context, score_type, metric_type) %>%
      dplyr::summarise(metric_sd = sd(metric_value,na.rm = T),
                       metric_value = mean(metric_value, na.rm = T)) %>%
      dplyr::ungroup() %>%
      tidyr::separate(col = context,
                      into = c("polarity","mz_class",
                               "analyzer", "source",
                               "Species", "Sample_type"),
                      sep = ",", remove = F)
  }
  
  if("cv_split" %in% colnames(annot_df)){
    colnames(annot_df)[colnames(annot_df) == "cv_split"] = "db"
  }
  
  if (!metrics_per_grp){
    annot_df = annot_df %>%
      dplyr::select(group_name, FDR_pct, msm_fdr_annots, pred_fdr_annots, db) %>%
      tidyr::separate(group_name,
                      into = c("ds_id", "adduct"), sep = ",",
                      remove = F) %>%
      dplyr::filter(ds_id %in% context_df$ds_id) %>% 
      dplyr::group_by(ds_id, FDR_pct,db) %>%
      dplyr::summarise(msm_fdr_per_ds = sum(msm_fdr_annots),
                       pred_fdr_per_ds = sum(pred_fdr_annots)) %>%
      dplyr::ungroup() %>%
      as.data.frame()
  }
  else{
    colnames(annot_df)[colnames(annot_df) == "msm_fdr_annots"] = "msm_fdr_per_ds"
    colnames(annot_df)[colnames(annot_df) == "pred_fdr_annots"] = "pred_fdr_per_ds"
  }
  annot_df = annot_df %>%
    dplyr::mutate(Difference = pred_fdr_per_ds - msm_fdr_per_ds,
                  diff_sign = ifelse(Difference < 0,-1,1),
                  LogDiff = diff_sign * log10(abs(Difference) + 1),
                  FC = (pred_fdr_per_ds + 1) / (msm_fdr_per_ds + 1),
                  LFC = log2(FC))
  if(!is.null(FDR.pct)){
    annot_df = annot_df %>%
      dplyr::filter(FDR_pct == FDR.pct)
  }
  if(metrics_per_grp){
    all_res = left_join(eval_df, annot_df, by = c("group_name")) %>% suppressWarnings()
  }
  else{
    if(mean_eval_context){
      eval_df$source_analyzer = paste0(eval_df$source, "_", eval_df$analyzer)
      return(eval_df)
    }
    all_res = left_join(eval_df, annot_df, by = c("ds_id")) %>% suppressWarnings()
  }
  all_res$source_analyzer = paste0(all_res$source, "_", all_res$analyzer)
  all_res$FDR_pct = as.character(all_res$FDR_pct)
  return(all_res)
}

prepare_enrich_table = function(enrich_res, ds_list,
                                sel_context,
                                filter_by_adjpval = T,
                                min_TP = 3, use_FE = T,
                                kingdom = c("Animal", "Plant")){
  
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
    dplyr::filter(!is.na(context)) %>% 
    dplyr::select(-context) %>%
    dplyr::distinct()
  
  return(plot_df)
  
}