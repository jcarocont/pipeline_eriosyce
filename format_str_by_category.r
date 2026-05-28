#!/usr/bin/env Rapp
#| description: Filtra un archivo STR por categoría de loci (RF, LFMM, OR, local-adapt, etc.)

#| description: Path al STR completo (tab-sep, con header)
str_file <- NULL

#| description: Path al RDS de clasificación de loci (output de selective_loci_lfmm-gf.r)
loci_rds <- NULL

#| description: Categoría a exportar: selec_rf | selec_lmmf | selected_or | selected_and | true_neutral | neutral_rf | neutral_lmmf
categoria <- NULL

#| description: Directorio de salida
outdir <- "."

library(dartR.base)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── validación ────────────────────────────────────────────────────────────────
if (is.null(str_file))  stop("--str_file es requerido")
if (is.null(loci_rds))  stop("--loci_rds es requerido")
if (is.null(categoria)) stop("--categoria es requerido")

loci_res <- readRDS(loci_rds)

categorias_validas <- names(loci_res)
if (!categoria %in% categorias_validas)
  stop(sprintf("Categoría '%s' no existe en el RDS.\nCategorías disponibles: %s",
               categoria, paste(categorias_validas, collapse = ", ")))

loci_vec <- loci_res[[categoria]]
cat(sprintf("Categoría: %s | Loci: %d\n", categoria, length(loci_vec)))

# ── cargar STR ────────────────────────────────────────────────────────────────
str_df <- read.csv(str_file, sep = "\t", check.names = FALSE)
cat(sprintf("STR cargado: %d filas × %d columnas\n", nrow(str_df), ncol(str_df)))

# ── normalizar nombres de columnas al formato del RDS (ej: 100264164-42-A/C) ──
col_loci <- colnames(str_df)[-1]
col_norm <- gsub("^X", "", col_loci)
col_norm <- gsub(
  pattern     = "(\\d+)\\.(\\d+)\\.(\\w)\\.(\\w)",
  replacement = "\\1-\\2-\\3/\\4",
  x           = col_norm
)

keep_idx <- which(col_norm %in% loci_vec)
cat(sprintf("Loci encontrados en STR: %d / %d\n", length(keep_idx), length(loci_vec)))

if (length(keep_idx) == 0)
  stop("Ningún locus de la categoría fue encontrado en el STR. Revisa el formato de nombres.")

# ── filtrar y exportar ────────────────────────────────────────────────────────
str_out <- str_df[, c(1, keep_idx + 1)]

outfile <- file.path(outdir, paste0("strfile_", categoria, ".str"))
write.table(str_out, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("✓ Guardado: %s (%d loci)\n", outfile, length(keep_idx)))
