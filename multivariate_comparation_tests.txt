#!/usr/bin/env Rapp
#| description: MRM, db-RDA y partición de varianza genómica

#| description: Path a la matriz de distancia genética (BayesAss, txt)
genetic_distance_file <- ""

#| description: Path al CSV de importancia de variables LFMM/LEA
varimportance_lea <- ""

#| description: Path al CSV de importancia de variables GF
varimportance_gf <- ""

#| description: Path al CSV de variables ambientales
env_file <- ""

#| description: Path al CSV de coordenadas (lon, lat)
coords_file <- ""

#| description: Path de salida para el RDS de resultados
outfile <- ""

#| description: Path de salida para el PNG del diagrama Euler
euler_png <- ""


  library(ecodist)
  library(geosphere)
  library(vegan)
  library(proxy)
  library(stringr)
  library(eulerr)
  library(grid)

# --- DISTANCIA GENÉTICA ---
genetic_distance <- read.csv(genetic_distance_file, header = FALSE, row.names = 1)
ronames <- rownames(genetic_distance)
genetic_distance <- genetic_distance[, -ncol(genetic_distance)]

splits <- str_split_fixed(ronames, "_", n = 3)[, 1:2]
unidos <- apply(splits, 1, paste, collapse = "_")
colnames(genetic_distance) <- unidos
rownames(genetic_distance) <- unidos

# --- VARIABLES AMBIENTALES ---
lea <- read.csv(varimportance_lea); lea$X <- NULL
gf  <- read.csv(varimportance_gf);  gf$X  <- NULL

lea <- lea[order(lea$n_sig_snps,        decreasing = TRUE), ]
gf  <- gf[order(gf$total_importance,    decreasing = TRUE), ]

selectlea <- lea[1:4, ]
selectgf  <- gf[1:4, ]

cat("\n=== Variables Seleccionadas ===\n")
cat("\nLFMM (LEA) - Top 4:\n");  print(selectlea[, c("variable", "n_sig_snps")])
cat("\nGF - Top 4:\n");           print(selectgf[,  c("variable", "total_importance")])

env_all  <- read.csv(env_file, row.names = 2)
env_all  <- env_all[unidos, ]
env_lfmm <- env_all[, selectlea$variable, drop = FALSE]
env_gf   <- env_all[, selectgf$variable,  drop = FALSE]
varnames <- unique(c(colnames(env_lfmm), colnames(env_gf)))
env_mix  <- env_all[, varnames, drop = FALSE]

cat("\n=== Variables Ambientales ===\n")
cat("LFMM:", ncol(env_lfmm), "variables\n")
cat("GF:",   ncol(env_gf),   "variables\n")
cat("Combinadas (únicas):", ncol(env_mix), "variables\n")

# --- DISTANCIA GEOGRÁFICA ---
coord_df <- read.csv(coords_file, row.names = 2)
coord_df <- coord_df[unidos, c("lon", "lat")]

geo_dist <- distm(coord_df, fun = distGeo) / 1000
rownames(geo_dist) <- unidos
colnames(geo_dist) <- unidos

# --- DISTANCIAS AMBIENTALES ---
env_lfmm_std <- scale(env_lfmm)
env_gf_std   <- scale(env_gf)
env_mix_std  <- scale(env_mix)

amb_lfmm_dist <- dist(env_lfmm_std, method = "euclidean")
amb_gf_dist   <- dist(env_gf_std,   method = "euclidean")
amb_mix_dist  <- dist(env_mix_std,  method = "euclidean")

# --- MANTEL ---
mantel_geo  <- mantel(as.dist(genetic_distance), as.dist(geo_dist),      permutations = 999)
mantel_lfmm <- mantel(as.dist(genetic_distance), as.dist(amb_lfmm_dist), permutations = 999)
mantel_gf   <- mantel(as.dist(genetic_distance), as.dist(amb_gf_dist),   permutations = 999)
mantel_mix  <- mantel(as.dist(genetic_distance), as.dist(amb_mix_dist),  permutations = 999)

cat("\nGen ~ Geo:\n");       cat("  r =", round(mantel_geo$statistic,  4), ", p =", mantel_geo$signif,  "\n")
cat("\nGen ~ Amb_LFMM:\n"); cat("  r =", round(mantel_lfmm$statistic, 4), ", p =", mantel_lfmm$signif, "\n")
cat("\nGen ~ Amb_GF:\n");   cat("  r =", round(mantel_gf$statistic,   4), ", p =", mantel_gf$signif,   "\n")
cat("\nGen ~ Amb_Mix:\n");  cat("  r =", round(mantel_mix$statistic,  4), ", p =", mantel_mix$signif,  "\n")

# --- MRM ---
gen_dist_vec <- as.vector(as.dist(genetic_distance))
geo_dist_vec <- as.vector(as.dist(geo_dist))
amb_lfmm_vec <- as.vector(amb_lfmm_dist)
amb_gf_vec   <- as.vector(amb_gf_dist)
amb_mix_vec  <- as.vector(amb_mix_dist)

mrm_lfmm        <- MRM(gen ~ geo + amb, data = data.frame(gen = gen_dist_vec, geo = geo_dist_vec, amb = amb_lfmm_vec), nperm = 999)
mrm_gf          <- MRM(gen ~ geo + amb, data = data.frame(gen = gen_dist_vec, geo = geo_dist_vec, amb = amb_gf_vec),   nperm = 999)
mrm_mix         <- MRM(gen ~ geo + amb, data = data.frame(gen = gen_dist_vec, geo = geo_dist_vec, amb = amb_mix_vec),  nperm = 999)
mrm_interaction <- MRM(gen ~ geo * amb, data = data.frame(gen = gen_dist_vec, geo = geo_dist_vec, amb = amb_mix_vec),  nperm = 999)

cat("\nMRM: Gen ~ Geo + Amb_LFMM\n"); print(mrm_lfmm)
cat("\nMRM: Gen ~ Geo + Amb_GF\n");   print(mrm_gf)
cat("\nMRM: Gen ~ Geo + Amb_Mix\n");  print(mrm_mix)
cat("\nMRM con interacción: Gen ~ Geo * Amb\n"); print(mrm_interaction)

# --- db-RDA ---
geo_pcoa <- cmdscale(geo_dist, k = min(10, nrow(coord_df) - 1), eig = TRUE)
geo_axes <- geo_pcoa$points
colnames(geo_axes) <- paste0("Geo", 1:ncol(geo_axes))

dbrda_lfmm <- capscale(genetic_distance ~ ., data = as.data.frame(env_lfmm_std))
dbrda_gf   <- capscale(genetic_distance ~ ., data = as.data.frame(env_gf_std))
dbrda_mix  <- capscale(genetic_distance ~ ., data = as.data.frame(env_mix_std))
dbrda_full <- capscale(genetic_distance ~ . + geo_axes, data = as.data.frame(env_mix_std))

anova_lfmm <- anova.cca(dbrda_lfmm, permutations = 999)
anova_gf   <- anova.cca(dbrda_gf,   permutations = 999)
anova_mix  <- anova.cca(dbrda_mix,  permutations = 999)
anova_full <- anova.cca(dbrda_full, permutations = 999)

cat("\ndb-RDA LFMM:\n");              cat("  R² ajustado:", round(RsquareAdj(dbrda_lfmm)$adj.r.squared, 4), "\n  p-value:", anova_lfmm$`Pr(>F)`[1], "\n")
cat("\ndb-RDA GF:\n");                cat("  R² ajustado:", round(RsquareAdj(dbrda_gf)$adj.r.squared,   4), "\n  p-value:", anova_gf$`Pr(>F)`[1],   "\n")
cat("\ndb-RDA Combinada (LFMM+GF):\n"); cat("  R² ajustado:", round(RsquareAdj(dbrda_mix)$adj.r.squared, 4), "\n  p-value:", anova_mix$`Pr(>F)`[1], "\n")
cat("\ndb-RDA Full:\n");              cat("  R² ajustado:", round(RsquareAdj(dbrda_full)$adj.r.squared, 4), "\n  p-value:", anova_full$`Pr(>F)`[1], "\n")

# --- db-RDA PARCIALES ---
cat("\n=== db-RDA Parciales ===\n")

dbrda_lfmm_partial <- capscale(genetic_distance ~ . + Condition(as.matrix(env_gf_std))   + Condition(geo_axes), data = as.data.frame(env_lfmm_std))
dbrda_gf_partial   <- capscale(genetic_distance ~ . + Condition(as.matrix(env_lfmm_std)) + Condition(geo_axes), data = as.data.frame(env_gf_std))
dbrda_geo_partial  <- capscale(genetic_distance ~ geo_axes + Condition(as.matrix(cbind(env_lfmm_std, env_gf_std))))

anova_lfmm_partial <- anova.cca(dbrda_lfmm_partial, permutations = 999)
anova_gf_partial   <- anova.cca(dbrda_gf_partial,   permutations = 999)
anova_geo_partial  <- anova.cca(dbrda_geo_partial,  permutations = 999)

cat("\nLFMM | GF+Geo:\n"); cat("  R² ajustado:", round(RsquareAdj(dbrda_lfmm_partial)$adj.r.squared, 4), "\n  p-value:", anova_lfmm_partial$`Pr(>F)`[1], "\n")
cat("\nGF | LFMM+Geo:\n"); cat("  R² ajustado:", round(RsquareAdj(dbrda_gf_partial)$adj.r.squared,   4), "\n  p-value:", anova_gf_partial$`Pr(>F)`[1],   "\n")
cat("\nGeo | LFMM+GF:\n"); cat("  R² ajustado:", round(RsquareAdj(dbrda_geo_partial)$adj.r.squared,  4), "\n  p-value:", anova_geo_partial$`Pr(>F)`[1],  "\n")

# --- VARPART ---
cat("\n=== Partición de Varianza (varpart) ===\n")
varpart_result <- varpart(genetic_distance, env_lfmm_std, env_gf_std, geo_axes)
print(varpart_result)

frac     <- varpart_result$part$indfract$Adj.R.square
frac_a   <- max(0, frac[1])
frac_b   <- max(0, frac[2])
frac_c   <- max(0, frac[3])
frac_ab  <- max(0, frac[4])
frac_ac  <- max(0, frac[5])
frac_bc  <- max(0, frac[6])
frac_abc <- max(0, frac[7])

cat("\n=== Fracciones de Varianza ===\n")
cat("  [a] LFMM único:",  round(frac_a   * 100, 2), "%\n")
cat("  [b] GF único:",    round(frac_b   * 100, 2), "%\n")
cat("  [c] Geo único:",   round(frac_c   * 100, 2), "%\n")
cat("  [ab] LFMM:GF:",    round(frac_ab  * 100, 2), "%\n")
cat("  [ac] LFMM:Geo:",   round(frac_ac  * 100, 2), "%\n")
cat("  [bc] GF:Geo:",     round(frac_bc  * 100, 2), "%\n")
cat("  [abc] Triple:",    round(frac_abc * 100, 2), "%\n")

# --- EULER ---
area_lfmm    <- max(0, (frac_a + frac_ab + frac_ac + frac_abc) * 100)
area_gf      <- max(0, (frac_b + frac_ab + frac_bc + frac_abc) * 100)
area_geo     <- max(0, (frac_c + frac_ac + frac_bc + frac_abc) * 100)
area_lfmm_gf  <- max(0, (frac_ab + frac_abc) * 100)
area_lfmm_geo <- max(0, (frac_ac + frac_abc) * 100)
area_gf_geo   <- max(0, (frac_bc + frac_abc) * 100)
area_triple   <- max(0, frac_abc * 100)

pct_explained   <- area_lfmm + area_gf + area_geo - area_lfmm_gf - area_lfmm_geo - area_gf_geo + area_triple
pct_unexplained <- max(0, 100 - pct_explained)

cat("\n=== Composición de Varianza (%) ===\n")
cat("  Explicada total:", round(pct_explained,   2), "%\n")
cat("  Residual:",        round(pct_unexplained, 2), "%\n")

combinations <- c(
  "LFMM"               = frac_a   * 100,
  "GF"                 = frac_b   * 100,
  "Geografía"          = frac_c   * 100,
  "LFMM&GF"            = frac_ab  * 100,
  "LFMM&Geografía"     = frac_ac  * 100,
  "GF&Geografía"       = frac_bc  * 100,
  "LFMM&GF&Geografía"  = frac_abc * 100
)

fit <- euler(combinations)
cat("\nError de ajuste Euler:", round(fit$stress * 100, 2), "%\n")

png(euler_png, width = 3000, height = 2400, res = 300)
plot(
  fit,
  fills     = list(fill = c("#87CEEB", "#FFB6C1", "#90EE90"), alpha = 0.6),
  edges     = list(col  = c("skyblue4", "pink4", "darkgreen"), lwd = 2),
  labels    = list(cex  = 1.2, font = 2),
  quantities = list(cex = 1.1, font = 1)
)
grid.text("Partición de Varianza Genómica (db-RDA)", y = unit(0.96, "npc"), gp = gpar(fontsize = 18, fontface = "bold"))
grid.text(
  sprintf("LFMM: %.2f%%\nGF: %.2f%%\nGeografía: %.2f%%\nResidual: %.2f%%", area_lfmm, area_gf, area_geo, pct_unexplained),
  x = unit(0.02, "npc"), y = unit(0.05, "npc"), just = "left", gp = gpar(fontsize = 11, col = "gray30")
)
dev.off()
cat(sprintf("\n✓ Diagrama guardado: %s\n", euler_png))

# --- GUARDAR RDS ---
results_list <- list(
  genetic_distance    = genetic_distance,
  geo_dist            = geo_dist,
  amb_lfmm_dist       = amb_lfmm_dist,
  amb_gf_dist         = amb_gf_dist,
  amb_mix_dist        = amb_mix_dist,
  env_lfmm_std        = env_lfmm_std,
  env_gf_std          = env_gf_std,
  env_mix_std         = env_mix_std,
  geo_axes            = geo_axes,
  mantel_geo          = mantel_geo,
  mantel_lfmm         = mantel_lfmm,
  mantel_gf           = mantel_gf,
  mantel_mix          = mantel_mix,
  mrm_lfmm            = mrm_lfmm,
  mrm_gf              = mrm_gf,
  mrm_mix             = mrm_mix,
  mrm_interaction     = mrm_interaction,
  dbrda_lfmm          = dbrda_lfmm,
  dbrda_gf            = dbrda_gf,
  dbrda_mix           = dbrda_mix,
  dbrda_full          = dbrda_full,
  dbrda_lfmm_partial  = dbrda_lfmm_partial,
  dbrda_gf_partial    = dbrda_gf_partial,
  dbrda_geo_partial   = dbrda_geo_partial,
  anova_lfmm          = anova_lfmm,
  anova_gf            = anova_gf,
  anova_mix           = anova_mix,
  anova_full          = anova_full,
  anova_lfmm_partial  = anova_lfmm_partial,
  anova_gf_partial    = anova_gf_partial,
  anova_geo_partial   = anova_geo_partial,
  varpart_result      = varpart_result,
  variance_fractions  = list(
    frac_a = frac_a, frac_b = frac_b, frac_c = frac_c,
    frac_ab = frac_ab, frac_ac = frac_ac, frac_bc = frac_bc, frac_abc = frac_abc,
    pct_explained = pct_explained, pct_unexplained = pct_unexplained
  ),
  euler_fit = fit
)

saveRDS(results_list, file = outfile)
cat(sprintf("✅ Resultados guardados: %s\n", outfile))
