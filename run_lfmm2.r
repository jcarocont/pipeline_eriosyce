#!/usr/bin/env Rapp
#| description: Asociación genotipo-ambiente con lfmm2. Clasifica SNPs en neutrales/selectivos y genera plots.

#| description: Ruta al archivo STR (tab-separated, sin header)
str_file <- NULL

#| description: Ruta al CSV de variables ambientales (con columnas id, ...)
env_file <- NULL

#| description: Número de factores latentes para lfmm2
k_latent <- 3L

#| description: Umbral de p-valor para clasificar SNPs como asociados
pval_threshold <- 0.05

#| description: Ruta al archivo RDS de salida del objeto lfmm2.test
output <- "asoc_ambiental.rds"

#| description: Directorio de salida para plots y tablas
outdir <- "."

library(LEA)
library(tidyverse)
library(ggplot2)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# 1. CARGAR Y PREPARAR STR
# ══════════════════════════════════════════════════════════════════════════════
cat("── Cargando STR:", str_file, "\n")
str_raw <- read.csv(str_file, sep = "\t", header = FALSE)

names_vector <- as.vector(na.omit(str_raw[1, 1:(ncol(str_raw) - 1)]))
rownames_vec <- str_raw[2:nrow(str_raw), 1]
str_df       <- str_raw[2:nrow(str_raw), 2:ncol(str_raw)]
colnames(str_df) <- names_vector

# una fila por individuo (haploid: cada 2 filas en STR)
ind_names <- unique(rownames_vec)
str_df1   <- str_df[seq(1, by = 2, length.out = length(ind_names)), ]
rownames(str_df1) <- ind_names
str_df1[str_df1 == -9] <- 0
str_df1[] <- lapply(str_df1, function(x) if (is.character(x)) as.numeric(x) else x)
str_mat <- as.matrix(str_df1)

cat("   Individuos:", nrow(str_mat), "| SNPs:", ncol(str_mat), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# 2. CARGAR Y PREPARAR VARIABLES AMBIENTALES
# ══════════════════════════════════════════════════════════════════════════════
cat("── Cargando variables ambientales:", env_file, "\n")
env_raw <- read.csv(env_file, sep = ",")
env_raw <- env_raw[env_raw$id %in% ind_names, ]
rownames(env_raw) <- env_raw$id
env_raw  <- env_raw[ind_names, ]          # mismo orden que STR
env_mat  <- as.matrix(env_raw[, 3:ncol(env_raw)])

cat("   Variables:", ncol(env_mat), "| Individuos:", nrow(env_mat), "\n")

# escalar
env_esc <- apply(env_mat, 2, function(x) (x - mean(x)) / sd(x))
cat("   Variables tras escalado:", colnames(env_esc), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# 3. LFMM2
# ══════════════════════════════════════════════════════════════════════════════
cat("── Ajustando lfmm2 con K =", k_latent, "factores latentes...\n")
latfact <- lfmm2(input = str_mat, env = env_esc, K = k_latent)

cat("── Corriendo lfmm2.test...\n")
assoc <- lfmm2.test(object = latfact, env = env_mat, input = str_mat)

cat("   GIF por variable:\n")
print(round(assoc$gif, 3))

saveRDS(assoc, file = output)
cat("   RDS guardado en:", output, "\n")

# ══════════════════════════════════════════════════════════════════════════════
# 4. CLASIFICAR SNPs
# ══════════════════════════════════════════════════════════════════════════════
cat("── Clasificando SNPs (umbral p <=", pval_threshold, ")...\n")

p_adj      <- assoc$pvalues
pvals_bool <- p_adj < pval_threshold

# limpiar nombres de columnas (quitar "Response X.")
valnames <- gsub("Response ", "", colnames(p_adj))
valnames <- unique(sub("\\..*$", "", sort(valnames)))
colnames(pvals_bool) <- valnames

# proporción de individuos asociados por SNP
col_any <- apply(pvals_bool, 2, mean)

neutral        <- names(col_any[col_any == 0])
some_selective <- names(col_any[col_any > 0])

cat("   SNPs neutrales:", length(neutral), "\n")
cat("   SNPs con alguna asociación:", length(some_selective), "\n")

writeLines(neutral,        file.path(outdir, "sitios_neutrales.txt"))
writeLines(some_selective, file.path(outdir, "sitios_selectivos.txt"))
cat("   Listas guardadas en:", outdir, "\n")

# ══════════════════════════════════════════════════════════════════════════════
# 5. TABLA DE IMPORTANCIA DE VARIABLES
# ══════════════════════════════════════════════════════════════════════════════
cat("── Calculando importancia de variables ambientales...\n")

sig_counts <- rowSums(p_adj < pval_threshold, na.rm = TRUE)
mean_abs_z <- rowMeans(abs(assoc$zscores), na.rm = TRUE)

var_importance <- tibble(
  variable   = rownames(p_adj),
  n_sig_snps = sig_counts,
  mean_abs_z = mean_abs_z,
  gif        = assoc$gif
) |>
  arrange(desc(n_sig_snps), desc(mean_abs_z))

print(var_importance)

write.csv(var_importance, file.path(outdir, "varimportance_lea.csv"), row.names = FALSE)
cat("   Tabla guardada en:", file.path(outdir, "varimportance_lea.csv"), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# 6. PLOTS GGPLOT
# ══════════════════════════════════════════════════════════════════════════════
cat("── Generando plots...\n")

# ── 6a. Histograma de proporción de asociación por SNP ───────────────────────
col_any_df <- tibble(prop_assoc = col_any[col_any > 0])

p_hist <- ggplot(col_any_df, aes(x = prop_assoc)) +
  geom_histogram(bins = 10, fill = "#4daf4a", color = "white", linewidth = 0.3) +
  labs(
    x     = "Proporción de individuos asociados",
    y     = "Número de SNPs",
    title = "Distribución de asociación SNP-ambiente",
    subtitle = paste0(
      length(some_selective), " SNPs con alguna asociación  |  ",
      length(neutral), " neutrales"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(outdir, "hist_prop_asociacion.png"),
  p_hist, width = 7, height = 5, dpi = 300
)
cat("   Guardado: hist_prop_asociacion.png\n")

# ── 6b. Importancia de variables (barplot) ────────────────────────────────────
p_var <- ggplot(var_importance, aes(x = reorder(variable, n_sig_snps), y = n_sig_snps)) +
  geom_col(fill = "#377eb8", color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = n_sig_snps),
    hjust = -0.2, size = 3.5
  ) +
  coord_flip() +
  labs(
    x     = "Variable ambiental",
    y     = paste0("N° SNPs significativos (p < ", pval_threshold, ")"),
    title = "Importancia de variables ambientales (lfmm2)"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(outdir, "importancia_variables.png"),
  p_var, width = 8, height = 5, dpi = 300
)
cat("   Guardado: importancia_variables.png\n")

# ── 6c. GIF por variable ──────────────────────────────────────────────────────
gif_df <- tibble(
  variable = names(assoc$gif),
  gif      = assoc$gif
)

p_gif <- ggplot(gif_df, aes(x = reorder(variable, gif), y = gif)) +
  geom_col(
    aes(fill = gif > 0.8 & gif < 1.2),
    color = "white", linewidth = 0.2, show.legend = FALSE
  ) +
  geom_hline(yintercept = c(0.8, 1.2), linetype = "dashed", color = "red", linewidth = 0.6) +
  scale_fill_manual(values = c("TRUE" = "#4daf4a", "FALSE" = "#e41a1c")) +
  coord_flip() +
  labs(
    x        = "Variable ambiental",
    y        = "GIF (Genomic Inflation Factor)",
    title    = "GIF por variable ambiental",
    subtitle = "Verde = rango aceptable (0.8 – 1.2)"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(outdir, "gif_variables.png"),
  p_gif, width = 8, height = 5, dpi = 300
)
cat("   Guardado: gif_variables.png\n")

cat("\n✓ Análisis completo.\n")
