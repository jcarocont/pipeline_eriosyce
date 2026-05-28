#!/home/agente/.local/bin/Rapp
# =============================================================================
# composite_turnover_full.R
# Análisis completo de turnover genómico + visualizaciones diagnósticas
# Integra: composite_turnover_merged.R + plots_gf.ipynb + curvas_gf.ipynb
# =============================================================================

library(vcfR)
library(adegenet)
library(ape)
library(LEA)
library(tidyverse)
library(furrr)
library(ranger)
library(rlang)
source("turnover_functions.R")

# ─────────────────────────────────────────────
# CONFIGURACIÓN
# ─────────────────────────────────────────────
NUM_BINS       <- 101
THR_FRAC_GOOD  <- 0.2   # turnover bien distribuido
THR_FRAC_BAD   <- 0.4   # turnover hiper-concentrado (posible sobreajuste)
THR_COR_GOOD   <- 0.6   # alta estabilidad de la curva
THR_COR_BAD    <- 0.4   # baja estabilidad de la curva

plan(multisession, workers = max(1, parallel::detectCores() - 1))
message(sprintf("[CONFIG] Núcleos disponibles: %d | NUM_BINS: %d",
                parallel::detectCores(), NUM_BINS))

# ─────────────────────────────────────────────
# 0. CARGA DE DATOS
# ─────────────────────────────────────────────
message("\n[0] Cargando datos y modelos...")

env <- read.csv("envar_final-corto.csv") %>%
  tibble::column_to_rownames("id") %>%
  select(-2, -X) %>%
  as.data.frame()

message(sprintf("[0] Tabla ambiental cargada: %d localidades × %d variables",
                nrow(env), ncol(env)))

model_list <- readRDS("forest_model.rds")
message(sprintf("[0] Modelos cargados: %d SNPs", length(model_list)))

# ─────────────────────────────────────────────
# 1. FILTRADO Y CACHE DE treeInfo
# ─────────────────────────────────────────────
message("\n[1] Filtrando modelos (r² > 0) y cacheando treeInfo...")

clean_models <- model_list %>%
  keep(~ inherits(.$model, "ranger") && .$model$r.squared > 0) %>%
  imap(~ { .x$model$.__treeInfo_cache__ <- treeInfo(.x$model); .x })

selective_sites <- names(clean_models)
neutral_sites   <- setdiff(names(model_list), selective_sites)

write.csv(data.frame(snp = selective_sites), "sitios_selectivos_rf.csv", row.names = FALSE)
write.csv(data.frame(snp = neutral_sites),   "sitios_neutros_rf.csv",    row.names = FALSE)

message(sprintf("[1] Sitios selectivos: %d | Sitios neutros: %d",
                length(selective_sites), length(neutral_sites)))

rm(model_list); gc()

# ─────────────────────────────────────────────
# 2A. RANKING GLOBAL DE VARIABLES
# ─────────────────────────────────────────────
message("\n[2A] Calculando ranking global de variables...")

compute_global_ranking <- function(models) {
  map_dfr(models, ~ {
    vi <- .x$model$variable.importance
    tibble(variable = names(vi), importance = as.numeric(vi))
  }, .id = "SNP_ID") %>%
    group_by(variable) %>%
    summarise(
      mean_imp   = mean(importance,   na.rm = TRUE),
      median_imp = median(importance, na.rm = TRUE),
      sd_imp     = sd(importance,     na.rm = TRUE),
      n_models   = n(),
      .groups    = "drop"
    ) %>%
    arrange(desc(mean_imp)) %>%
    mutate(rank = row_number())
}

ranking_global <- compute_global_ranking(clean_models)
message(sprintf("[2A] Variables rankeadas: %d", nrow(ranking_global)))
print(ranking_global)
write.csv(ranking_global, "ranking_global_variables.csv", row.names = FALSE)

# ─────────────────────────────────────────────
# 2B. TOP-3 VARIABLES POR SNP + PLOTS
# ─────────────────────────────────────────────
message("\n[2B] Calculando ranking top-3 por SNP...")

compute_top3_ranking <- function(models, n = 3) {
  map_dfr(models, ~ {
    vi <- sort(.x$model$variable.importance, decreasing = TRUE)
    tibble(
      var     = names(vi)[1:n],
      metrica = as.numeric(vi)[1:n],
      rank    = 1:n
    )
  }, .id = "snp_name") %>%
    mutate(metrica = as.numeric(metrica))
}

df_top3 <- compute_top3_ranking(clean_models)
message(sprintf("[2B] Filas generadas en top-3: %d", nrow(df_top3)))

p_count_all <- df_top3 %>%
  count(var) %>%
  arrange(n) %>%
  ggplot(aes(fct_reorder(var, n), n)) +
  geom_col() + coord_flip() +
  labs(title = "Frecuencia en top-3 (todas)", x = "Variable", y = "n") +
  theme_minimal()

p_count_filtered <- df_top3 %>%
  filter(metrica >= 1) %>%
  count(var) %>%
  arrange(n) %>%
  ggplot(aes(fct_reorder(var, n), n)) +
  geom_col() + coord_flip() +
  labs(title = "Frecuencia en top-3 (metrica >= 1)", x = "Variable", y = "n") +
  theme_minimal()

p_weighted_all <- df_top3 %>%
  group_by(var) %>%
  summarise(n = sum(metrica, na.rm = TRUE), .groups = "drop") %>%
  arrange(n) %>%
  ggplot(aes(fct_reorder(var, n), n)) +
  geom_col() + coord_flip() +
  labs(title = "Importancia ponderada top-3 (todas)", x = "Variable", y = "Suma importancia") +
  theme_minimal()

p_weighted_filtered <- df_top3 %>%
  filter(metrica >= 1) %>%
  group_by(var) %>%
  summarise(n = sum(metrica, na.rm = TRUE), .groups = "drop") %>%
  arrange(n) %>%
  ggplot(aes(fct_reorder(var, n), n)) +
  geom_col() + coord_flip() +
  labs(title = "Importancia ponderada top-3 (metrica >= 1)", x = "Variable", y = "Suma importancia") +
  theme_minimal()

message("[2B] Imprimiendo plots de top-3...")
print(p_count_all)
print(p_count_filtered)
print(p_weighted_all)
print(p_weighted_filtered)

# ─────────────────────────────────────────────
# 3. FUNCIONES DE TURNOVER
# ─────────────────────────────────────────────
message("\n[3] Definiendo funciones de turnover...")

extract_r2_contribution <- function(model, predictor_name) {
  r2 <- model$r.squared
  if (is.null(r2) || r2 <= 0) return(0)
  vi <- model$variable.importance
  if (is.null(vi) || !(predictor_name %in% names(vi))) return(0)
  r2 * (vi[predictor_name] / sum(vi, na.rm = TRUE))
}

extract_split_density_cached <- function(treeInfo_cache, x_vec, predictor_name, num_bins = NUM_BINS) {
  splitvals <- treeInfo_cache %>%
    filter(splitvarName == predictor_name) %>%
    pull(splitval)
  if (length(splitvals) < 2) return(NULL)
  rng    <- range(x_vec, na.rm = TRUE)
  grid_x <- seq(rng[1], rng[2], length.out = num_bins)
  tibble(
    x           = grid_x,
    dens_splits = density(splitvals, from = rng[1], to = rng[2], n = num_bins)$y,
    dens_data   = density(x_vec,     from = rng[1], to = rng[2], n = num_bins)$y
  )
}

calculate_turnover_ranger <- function(snp_entry, env_df, predictor_name, num_bins = NUM_BINS) {
  model <- snp_entry$model
  ti    <- model$.__treeInfo_cache__
  if (is.null(ti)) ti <- treeInfo(model)
  x_vec <- env_df[[predictor_name]]
  dens  <- extract_split_density_cached(ti, x_vec, predictor_name, num_bins)
  if (is.null(dens)) return(tibble(x = NA_real_, F_x = NA_real_)[0, ])
  r2_c  <- extract_r2_contribution(model, predictor_name)
  raw   <- cumsum(dens$dens_splits)
  F_x   <- raw / max(raw, na.rm = TRUE) * r2_c
  tibble(x = dens$x, F_x = F_x)
}

compute_all_turnover <- function(models, env_df, predictor_name, num_bins = NUM_BINS) {
  map_dfr(names(models), function(id) {
    snp_entry        <- models[[id]]
    snp_entry$SNP_ID <- id
    calculate_turnover_ranger(snp_entry, env_df, predictor_name, num_bins)
  }, .id = "SNP_ID", .progress = TRUE)
}

message("[3] Funciones de turnover definidas OK")

# ─────────────────────────────────────────────
# 4. CÁLCULO DE CURVAS DE TURNOVER
# ─────────────────────────────────────────────
target_vars <- ranking_global$variable[4:10]
message(sprintf("\n[4] Calculando turnover para %d variables: %s",
                length(target_vars), paste(target_vars, collapse = ", ")))

all_turnovers <- map(target_vars, ~ {
  message(sprintf("  [4] Procesando variable: %s", .x))
  result <- compute_all_turnover(clean_models, env, .x, NUM_BINS)
  message(sprintf("  [4] Done: %s — filas generadas: %d | SNPs con datos: %d",
                  .x, nrow(result), n_distinct(result$SNP_ID)))
  result
}) %>% set_names(target_vars)

message(sprintf("[4] Turnover calculado para %d variables.", length(all_turnovers)))

# ─────────────────────────────────────────────
# 5. EXPORTACIÓN DE CURVAS (CSV)
# ─────────────────────────────────────────────
message("\n[5] Exportando curvas de turnover a CSV...")

walk2(all_turnovers, names(all_turnovers), ~ {
  if (nrow(.x) > 0) {
    fname <- paste0("turnover_curves_", .y, ".csv")
    write.csv(.x, fname, row.names = FALSE)
    message(sprintf("  [5] Exportado: %s (%d filas)", fname, nrow(.x)))
  } else {
    message(sprintf("  [5] WARN: sin datos para %s, CSV omitido.", .y))
  }
})

# ─────────────────────────────────────────────
# 6. CONSTRUCCIÓN DE turnover_full
#    (equivalente al merge de los notebooks)
# ─────────────────────────────────────────────
message("\n[6] Construyendo turnover_full desde all_turnovers...")

turnover_full <- bind_rows(
  imap(all_turnovers, ~ mutate(.x, variable = .y))
) %>%
  left_join(
    map_dfr(clean_models, ~ {
      vi <- .x$model$variable.importance
      tibble(variable = names(vi), importance = as.numeric(vi),
             r_squared = .x$model$r.squared)
    }, .id = "SNP_ID"),
    by = c("SNP_ID", "variable")
  )

message(sprintf("[6] turnover_full: %d filas | %d SNPs | %d variables",
                nrow(turnover_full),
                n_distinct(turnover_full$SNP_ID),
                n_distinct(turnover_full$variable)))

# ─────────────────────────────────────────────
# 7. TABLAS RESUMEN
# ─────────────────────────────────────────────
message("\n[7] Calculando tablas resumen...")

# 7a. Resumen por SNP
snp_summary <- turnover_full |>
  group_by(SNP_ID) |>
  summarise(
    total_turnover = sum(F_x, na.rm = TRUE),
    mean_r2        = mean(r_squared, na.rm = TRUE),
    .groups        = "drop"
  )
message(sprintf("[7a] snp_summary: %d filas", nrow(snp_summary)))

# 7b. Matriz SNP × Variable
var_snp_matrix <- turnover_full |>
  group_by(SNP_ID, variable) |>
  summarise(total_turnover = sum(F_x, na.rm = TRUE), .groups = "drop")
message(sprintf("[7b] var_snp_matrix: %d filas", nrow(var_snp_matrix)))

# 7c. Importancia por variable
var_importance_summary <- turnover_full |>
  distinct(SNP_ID, variable, importance) |>
  group_by(variable) |>
  summarise(
    total_importance = sum(importance, na.rm = TRUE),
    n_snps           = n(),
    .groups          = "drop"
  ) |>
  mutate(percent_snps = 100 * n_snps / n_distinct(turnover_full$SNP_ID))
message(sprintf("[7c] var_importance_summary: %d variables", nrow(var_importance_summary)))

# 7d. Top-4 variables más importantes
top4_vars <- var_importance_summary |>
  slice_max(total_importance, n = 4) |>
  pull(variable)
message(sprintf("[7d] Top-4 variables: %s", paste(top4_vars, collapse = ", ")))

# 7e. Concentración del turnover (diagnóstico sobreajuste)
turnover_concentration <- turnover_full |>
  group_by(SNP_ID, variable) |>
  summarise(
    total_fx  = sum(F_x),
    peak_fx   = max(F_x),
    frac_peak = peak_fx / total_fx,
    .groups   = "drop"
  )
message(sprintf("[7e] turnover_concentration: %d filas", nrow(turnover_concentration)))

# 7f. Correlación de cada SNP con la curva media (estabilidad)
mean_curve <- turnover_full |>
  group_by(variable, x) |>
  summarise(mean_fx = mean(F_x, na.rm = TRUE), .groups = "drop")

curve_cor <- turnover_full |>
  left_join(mean_curve, by = c("variable", "x")) |>
  group_by(variable, SNP_ID) |>
  summarise(
    cor_to_mean = cor(F_x, mean_fx, use = "complete.obs"),
    .groups     = "drop"
  )
message(sprintf("[7f] curve_cor: %d filas", nrow(curve_cor)))

# 7g. Índice de zona (good / borderline / overfit) por SNP
curve_xy <- curve_cor |>
  left_join(turnover_concentration, by = c("SNP_ID", "variable")) |>
  group_by(SNP_ID) |>
  summarise(
    frac_peak   = max(frac_peak,   na.rm = TRUE),
    cor_to_mean = mean(cor_to_mean, na.rm = TRUE),
    .groups     = "drop"
  ) |>
  mutate(
    zone = case_when(
      frac_peak < THR_FRAC_GOOD & cor_to_mean > THR_COR_GOOD ~ "good",
      frac_peak > THR_FRAC_BAD  & cor_to_mean < THR_COR_BAD  ~ "overfit",
      TRUE ~ "borderline"
    )
  )

zone_counts <- curve_xy |>
  distinct(SNP_ID, zone) |>
  count(zone) |>
  mutate(pct = 100 * n / sum(n))

message("[7g] Distribución de zonas diagnósticas:")
print(zone_counts)

# SNPs a conservar (filtrando sobreajuste)
snps_keep <- curve_cor |>
  left_join(turnover_concentration, by = c("SNP_ID", "variable")) |>
  filter(frac_peak < 0.3, cor_to_mean > 0.5) |>
  distinct(SNP_ID)
message(sprintf("[7g] SNPs que pasan filtro diagnóstico (frac_peak<0.3 & cor>0.5): %d",
                nrow(snps_keep)))

# 7h. Exportar tablas resumen
write.csv(snp_summary,            "snp_summary.csv",            row.names = FALSE)
write.csv(var_importance_summary, "var_importance_summary.csv", row.names = FALSE)
write.csv(var_snp_matrix,         "var_snp_matrix.csv",         row.names = FALSE)
write.csv(curve_xy,               "snps_by_zone.csv",           row.names = FALSE)
write.csv(snps_keep,              "snps_diagnostico_ok.csv",    row.names = FALSE)
message("[7h] Tablas resumen exportadas.")

# ─────────────────────────────────────────────
# 8. PLOTS — CURVAS DE TURNOVER POR VARIABLE
#    (estilo composite_turnover_merged)
# ─────────────────────────────────────────────
message("\n[8] Generando plots de curvas de turnover por variable...")

plots_turnover <- imap(all_turnovers, ~ {
  df <- .x %>% filter(!is.na(x))
  if (nrow(df) == 0) {
    message(sprintf("  [8] WARN: sin datos para variable '%s', plot omitido.", .y))
    return(NULL)
  }
  gradient_curve <- df %>%
    group_by(x) %>%
    summarise(F_x_sum = sum(F_x, na.rm = TRUE), .groups = "drop")
  p <- ggplot() +
    geom_line(data = df,
              aes(x = x, y = F_x, group = SNP_ID),
              alpha = 0.12, color = "gray50") +
    geom_line(data = gradient_curve,
              aes(x = x, y = F_x_sum),
              linewidth = 1.2, color = "#0072B2") +
    labs(
      title    = paste("Turnover acumulado —", .y),
      subtitle = sprintf("N SNPs = %d", n_distinct(df$SNP_ID)),
      x        = .y,
      y        = expression(Sigma * R^2)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  ggsave(paste0("turnover_acumulado_", .y, ".png"), plot = p,
         width = 8, height = 6, dpi = 300)
  message(sprintf("  [8] Plot guardado: turnover_acumulado_%s.png", .y))
  p
})

plots_valid <- keep(plots_turnover, ~ !is.null(.x))
message(sprintf("[8] Plots de curvas: %d / %d generados.", length(plots_valid), length(target_vars)))
walk(plots_valid, print)

# ─────────────────────────────────────────────
# 9. PLOTS — IMPORTANCIA Y DISTRIBUCIÓN DE SNPs
# ─────────────────────────────────────────────
message("\n[9] Generando plots de importancia de variables...")

# 9a. SNPs con mayor recambio ambiental
p9a <- ggplot(snp_summary,
              aes(x = reorder(SNP_ID, total_turnover), y = total_turnover)) +
  geom_col() + coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_blank()) +
  labs(x = "SNP", y = "Turnover total", title = "SNPs con mayor recambio ambiental")
ggsave("snpxrecambio.png", plot = p9a, width = 8, height = 6, dpi = 300)
message("[9a] snpxrecambio.png guardado")
print(p9a)

# 9b. Heatmap SNP × Variable
p9b <- ggplot(var_snp_matrix,
              aes(x = variable, y = SNP_ID, fill = total_turnover)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Variable", y = "SNP", fill = "Turnover", title = "Turnover SNP × Variable")
ggsave("turnoverxsnp.png", plot = p9b, width = 8, height = 6, dpi = 300)
message("[9b] turnoverxsnp.png guardado")
print(p9b)

# 9c. Importancia total por variable
p9c <- ggplot(var_importance_summary,
              aes(x = reorder(variable, total_importance), y = total_importance)) +
  geom_col() + coord_flip() +
  theme_bw() +
  labs(x = "Variable ambiental", y = "Importancia total",
       title = "Importancia global de variables ambientales")
ggsave("importancia_global_varamb.png", plot = p9c, width = 8, height = 6, dpi = 300)
message("[9c] importancia_global_varamb.png guardado")
print(p9c)

# 9d. % de SNPs por variable
p9d <- ggplot(var_importance_summary,
              aes(x = reorder(variable, percent_snps), y = percent_snps)) +
  geom_col() + coord_flip() +
  theme_classic() +
  labs(x = "Variable", y = "% de SNPs", title = "Proporción de SNPs asociados por variable")
ggsave("importancia_snp_varamb.png", plot = p9d, width = 8, height = 6, dpi = 300)
message("[9d] importancia_snp_varamb.png guardado")
print(p9d)

# ─────────────────────────────────────────────
# 10. PLOTS — ENVELOPES Y QUANTILES
# ─────────────────────────────────────────────
message("\n[10] Generando plots de envelopes y cuantiles...")

# 10a. Envelope crudo (todas las curvas individuales)
p10a <- ggplot(turnover_full,
               aes(x = x, y = F_x, group = SNP_ID)) +
  geom_line(alpha = 0.02, linewidth = 0.2) +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  labs(x = "Gradiente ambiental", y = "Turnover",
       title = "Envelope de curvas de turnover por variable")
ggsave("envelope_turnover.png", plot = p10a, width = 10, height = 8, dpi = 300)
message("[10a] envelope_turnover.png guardado")
print(p10a)

# 10b. Envelope cuantílico (Q5/Q25/Q50/Q75/Q95)
quant_curve <- turnover_full |>
  group_by(variable, x) |>
  summarise(
    q05 = quantile(F_x, 0.05, na.rm = TRUE),
    q25 = quantile(F_x, 0.25, na.rm = TRUE),
    q50 = quantile(F_x, 0.50, na.rm = TRUE),
    q75 = quantile(F_x, 0.75, na.rm = TRUE),
    q95 = quantile(F_x, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

p10b <- ggplot(quant_curve, aes(x = x)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "grey80") +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "grey60") +
  geom_line(aes(y = q50), color = "black") +
  facet_wrap(~ variable, scales = "free") +
  theme_bw() +
  labs(x = "Gradiente ambiental", y = "Turnover",
       title = "Envelope cuantílico del turnover (diagnóstico de sobreajuste)")
ggsave("turnover_rangos_varianza.png", plot = p10b, width = 10, height = 8, dpi = 300)
message("[10b] turnover_rangos_varianza.png guardado")
print(p10b)

# ─────────────────────────────────────────────
# 11. PLOTS — TOP SNPs Y CURVAS GLOBALES
# ─────────────────────────────────────────────
message("\n[11] Generando plots de top-SNPs y curvas globales...")

# 11a. Curvas de los top-5 SNPs adaptativos
top_snps <- snp_summary |>
  slice_max(total_turnover, n = 5) |>
  pull(SNP_ID)
message(sprintf("[11a] Top-5 SNPs: %s", paste(top_snps, collapse = ", ")))

tf_top_snps <- turnover_full |> filter(SNP_ID %in% top_snps)

p11a <- ggplot(tf_top_snps,
               aes(x = x, y = F_x, color = variable)) +
  geom_line() +
  facet_wrap(~ SNP_ID, scales = "free_x") +
  theme_bw() +
  labs(x = "Gradiente", y = "Turnover",
       title = "Curvas de turnover por SNP (top adaptativos)")
ggsave("top_snps_turnover.png", plot = p11a, width = 10, height = 8, dpi = 300)
message("[11a] top_snps_turnover.png guardado")
print(p11a)

# 11b. Turnover promedio por variable
turnover_mean <- turnover_full |>
  group_by(variable, x) |>
  summarise(mean_turnover = mean(F_x, na.rm = TRUE), .groups = "drop")

p11b <- ggplot(turnover_mean,
               aes(x = x, y = mean_turnover)) +
  geom_line(color = "black", linewidth = 1.1) +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  labs(x = "Gradiente ambiental", y = "Turnover promedio por SNP",
       title = "Turnover promedio por variable")
ggsave("turnover_promedio.png", plot = p11b, width = 10, height = 8, dpi = 300)
message("[11b] turnover_promedio.png guardado")
print(p11b)

# 11c. Turnover acumulado global (suma) por variable
sum_curve <- turnover_full |>
  group_by(variable, x) |>
  summarise(sum_fx = sum(F_x, na.rm = TRUE), .groups = "drop")

p11c <- ggplot(sum_curve, aes(x = x, y = sum_fx)) +
  geom_line(color = "black", linewidth = 1.2) +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  labs(x = "Gradiente ambiental", y = "Turnover genómico acumulado",
       title = "Turnover genómico acumulado por variable")
ggsave("turnover_acumulado_global.png", plot = p11c, width = 10, height = 8, dpi = 300)
message("[11c] turnover_acumulado_global.png guardado")
print(p11c)

# 11d. Curvas individuales + global superpuestas (top-4 variables)
tf_top4 <- turnover_full |> filter(variable %in% top4_vars)

global_curve_top4 <- tf_top4 |>
  group_by(variable, x) |>
  summarise(F_x_global = sum(F_x, na.rm = TRUE), .groups = "drop")

p11d <- ggplot() +
  geom_line(data = tf_top4,
            aes(x = x, y = F_x, group = SNP_ID),
            color = "grey70", linewidth = 0.3, alpha = 0.4) +
  geom_line(data = global_curve_top4,
            aes(x = x, y = F_x_global),
            color = "black", linewidth = 1.2) +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  labs(x = "Gradiente ambiental", y = "Turnover genómico",
       title = "Turnover por SNP y respuesta genómica global (top-4 variables)")
ggsave("turnover_global_top4.png", plot = p11d, width = 10, height = 8, dpi = 300)
message("[11d] turnover_global_top4.png guardado")
print(p11d)

# ─────────────────────────────────────────────
# 12. PLOTS — DIAGNÓSTICO DE SOBREAJUSTE
# ─────────────────────────────────────────────
message("\n[12] Generando plots diagnósticos de sobreajuste...")

# 12a. Histograma de concentración (frac_peak)
p12a <- ggplot(turnover_concentration, aes(x = frac_peak)) +
  geom_histogram(bins = 40, fill = "grey30", color = "white") +
  theme_bw() +
  labs(x = "Fracción del turnover en el bin máximo", y = "Número de SNPs",
       title = "Concentración del turnover (diagnóstico de sobreajuste)")
ggsave("diagnostico_frac_peak_hist.png", plot = p12a, width = 8, height = 6, dpi = 300)
message("[12a] diagnostico_frac_peak_hist.png guardado")
print(p12a)

# 12b. Boxplot de frac_peak por variable
p12b <- ggplot(turnover_concentration,
               aes(x = variable, y = frac_peak)) +
  geom_boxplot(outlier.color = "firebrick") +
  coord_flip() +
  theme_bw() +
  labs(x = "Variable ambiental", y = "Frac. del turnover en el pico máximo",
       title = "Concentración del turnover por variable")
ggsave("diagnostico_frac_peak_boxplot.png", plot = p12b, width = 8, height = 6, dpi = 300)
message("[12b] diagnostico_frac_peak_boxplot.png guardado")
print(p12b)

# 12c. Correlación SNP ↔ curva media (estabilidad)
p12c <- ggplot(curve_cor,
               aes(x = variable, y = cor_to_mean)) +
  geom_boxplot(outlier.color = "firebrick") +
  coord_flip() +
  theme_bw() +
  labs(x = "Variable ambiental", y = "Correlación SNP ↔ curva media",
       title = "Estabilidad de las respuestas de SNPs (diagnóstico GF)")
ggsave("diagnostico_estabilidad_cor.png", plot = p12c, width = 8, height = 6, dpi = 300)
message("[12c] diagnostico_estabilidad_cor.png guardado")
print(p12c)

# 12d. Scatter diagnóstico: frac_peak vs cor_to_mean con zonas coloreadas
p12d <- curve_xy |>
  ggplot(aes(x = frac_peak, y = cor_to_mean)) +
  # Zona buena
  geom_rect(xmin = 0, xmax = THR_FRAC_GOOD,
            ymin = THR_COR_GOOD, ymax = 1,
            fill = "darkseagreen2", alpha = 0.3, inherit.aes = FALSE) +
  # Zona borderline
  geom_rect(xmin = THR_FRAC_GOOD, xmax = THR_FRAC_BAD,
            ymin = THR_COR_BAD, ymax = THR_COR_GOOD,
            fill = "khaki1", alpha = 0.3, inherit.aes = FALSE) +
  # Zona mala
  geom_rect(xmin = THR_FRAC_BAD, xmax = 1,
            ymin = 0, ymax = THR_COR_BAD,
            fill = "indianred1", alpha = 0.3, inherit.aes = FALSE) +
  geom_point(aes(color = zone), alpha = 0.8, size = 0.8) +
  scale_color_manual(values = c(good = "darkgreen", borderline = "goldenrod3", overfit = "firebrick")) +
  geom_vline(xintercept = c(THR_FRAC_GOOD, THR_FRAC_BAD), linetype = "dashed") +
  geom_hline(yintercept = c(THR_COR_GOOD, THR_COR_BAD), linetype = "dashed") +
  theme_bw() +
  labs(x = "Concentración del turnover (frac_peak)",
       y = "Correlación con curva media",
       title = "Diagnóstico de sobreajuste en Gradient Forest",
       subtitle = "Zonas indican estabilidad del ajuste SNP-específico",
       color = "Zona")
ggsave("diagnostico_scatter_zonas.png", plot = p12d, width = 8, height = 6, dpi = 300)
message("[12d] diagnostico_scatter_zonas.png guardado")
print(p12d)

# 12e. Dominancia de SNPs por variable
dominance_df <- turnover_full |>
  group_by(variable, SNP_ID) |>
  summarise(total_fx = sum(F_x), .groups = "drop") |>
  group_by(variable) |>
  mutate(frac = total_fx / sum(total_fx)) |>
  ungroup()

p12e <- dominance_df |>
  group_by(variable) |>
  summarise(top_frac = max(frac)) |>
  ggplot(aes(x = reorder(variable, top_frac), y = top_frac)) +
  geom_col(fill = "firebrick") +
  coord_flip() +
  theme_bw() +
  labs(x = "Variable ambiental",
       y = "Fracción del turnover explicada por el SNP dominante",
       title = "Dominancia de SNPs en el turnover (diagnóstico de sobreajuste)")
ggsave("diagnostico_dominancia_snp.png", plot = p12e, width = 8, height = 6, dpi = 300)
message("[12e] diagnostico_dominancia_snp.png guardado")
print(p12e)

# 12f. r² vs turnover total por SNP
p12f <- snp_summary |>
  ggplot(aes(x = mean_r2, y = total_turnover)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE, color = "#0072B2") +
  theme_bw() +
  labs(x = "R² promedio", y = "Turnover total",
       title = "Relación R² ↔ turnover total por SNP")
ggsave("r2_vs_turnover.png", plot = p12f, width = 8, height = 6, dpi = 300)
message("[12f] r2_vs_turnover.png guardado")
print(p12f)
snps_putativos <- snp_summary |>
  filter(
    mean_r2        > quantile(mean_r2,        0.50, na.rm = TRUE),
    total_turnover > quantile(total_turnover, 0.30, na.rm = TRUE),
    total_turnover < quantile(total_turnover, 0.70, na.rm = TRUE)
  )

message(sprintf("[putative] SNPs filtrados: %d", nrow(snps_putativos)))
write.csv(snps_putativos, "snps_putativos.csv", row.names = FALSE)
# ─────────────────────────────────────────────
# FIN
# ─────────────────────────────────────────────
message("\n========================================")
message("[FIN] Análisis completo terminado.")
message(sprintf("  SNPs selectivos analizados : %d", length(selective_sites)))
message(sprintf("  Variables de turnover      : %d (%s)", length(target_vars),
                paste(target_vars, collapse = ", ")))
message(sprintf("  SNPs que pasan diagnóstico : %d", nrow(snps_keep)))
message(sprintf("  Plots generados            : 12 bloques (PNG + pantalla)"))
message("  Archivos CSV exportados: snp_summary, var_importance_summary,")
message("    var_snp_matrix, snps_by_zone, snps_diagnostico_ok,")
message("    turnover_curves_<variable>, ranking_global_variables")
message("========================================\n")
