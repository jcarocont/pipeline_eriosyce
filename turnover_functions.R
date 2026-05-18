#!/usr/bin/env Rscript

# Funciones corregidas para trabajar con MATRICES

extract_r2_contribution <- function(model, predictor_name) {
  r_squared_total <- model$r.squared
  if (r_squared_total <= 0) return(0)
  
  var_importance <- model$variable.importance
  imp_pred <- var_importance[predictor_name]
  
  if (is.null(imp_pred) || is.na(imp_pred) || imp_pred <= 0) return(0)
  
  sum_importance <- sum(var_importance)
  r_squared_pred <- r_squared_total * (imp_pred / sum_importance)
  
  return(r_squared_pred)
}


extract_split_density <- function(model, predictor_name, env_matrix, num_bins = 101) {
  
  # Extraer splits
  all_splits <- lapply(1:model$num.trees, function(i) {
    tree <- ranger::treeInfo(model, tree = i)
    tree %>% 
      dplyr::filter(splitvarName == predictor_name) %>% 
      dplyr::select(splitval)
  }) %>% dplyr::bind_rows()
  
  if(nrow(all_splits) == 0) return(NULL)
  
  # *** CORRECCIÓN: Acceso correcto a columnas de matriz ***
  if (!predictor_name %in% colnames(env_matrix)) {
    warning(paste("Variable", predictor_name, "no existe en env"))
    return(NULL)
  }
  
  raw_data_values <- env_matrix[, predictor_name]
  range_val <- range(raw_data_values, na.rm = TRUE)
  grid_x <- seq(range_val[1], range_val[2], length.out = num_bins)
  
  # Densidad de splits
  dens_splits <- density(all_splits$splitval, 
                         from = range_val[1], 
                         to = range_val[2], 
                         n = num_bins)$y
  
  # Densidad de datos
  dens_data <- density(raw_data_values, 
                       from = range_val[1], 
                       to = range_val[2], 
                       n = num_bins)$y
  
  return(list(grid_x = grid_x, dens_splits = dens_splits, dens_data = dens_data))
}


calculate_turnover_rate <- function(density_list) {
  epsilon <- 1e-6
  turnover_rate <- density_list$dens_splits / (density_list$dens_data + epsilon)
  return(turnover_rate)
}


normalize_cumulative_turnover <- function(turnover_rate, r_squared_pred) {
  cumulative_turnover <- cumsum(turnover_rate)
  
  if (max(cumulative_turnover) == 0 || r_squared_pred == 0) {
    return(rep(0, length(turnover_rate)))
  }
  
  scale_factor <- r_squared_pred / max(cumulative_turnover)
  F_x <- cumulative_turnover * scale_factor
  
  return(F_x)
}


# *** FUNCIÓN PRINCIPAL CORREGIDA ***
calculate_turnover_ranger <- function(model, env, var, n_bins = 101) {
  
  r_squared_pred <- extract_r2_contribution(model, var)
  
  if (r_squared_pred == 0) {
    return(data.frame(x = NA, F_x = NA, R2_contribution = 0))
  }
  
  density_list <- extract_split_density(model, var, env, n_bins)
  
  if (is.null(density_list)) {
    warning(paste("Predictor", var, "sin splits"))
    return(data.frame(x = NA, F_x = NA, R2_contribution = r_squared_pred))
  }
  
  turnover_rate <- calculate_turnover_rate(density_list)
  F_x <- normalize_cumulative_turnover(turnover_rate, r_squared_pred)
  
  return(data.frame(x = density_list$grid_x, F_x = F_x, R2_contribution = r_squared_pred))
}
