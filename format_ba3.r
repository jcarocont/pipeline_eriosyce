#!/usr/bin/env Rapp
#| description: Convierte un STR filtrado por categoría al formato de input de BayesAss3

#| description: Path al STR (tab-separated, con header)
str_file <- NULL

#| description: Path al CSV de tabla de poblaciones (columnas: genotipo, popname_str)
poptable_file <- NULL

#| description: Path de salida para el archivo BayesAss3
outfile <- "bayesass3_input.txt"

library(tidyverse)

if (is.null(str_file))    stop("--str_file es requerido")
if (is.null(poptable_file)) stop("--poptable_file es requerido")

# ── cargar STR ────────────────────────────────────────────────────────────────
cat("── Cargando STR:", str_file, "\n")
df <- read.csv(str_file, sep = "\t")
cat(sprintf("   %d filas × %d columnas\n", nrow(df), ncol(df)))

# ── pivot a formato largo ─────────────────────────────────────────────────────
df_long <- df |>
  pivot_longer(cols = -genotipo, names_to = "locID", values_to = "coding") |>
  unique()

# ── cargar poptable ───────────────────────────────────────────────────────────
cat("── Cargando poptable:", poptable_file, "\n")
pops <- read.csv(poptable_file)
if ("X" %in% colnames(pops)) pops$X <- NULL
pops <- pops[, c("pop","id")]
pops <- pops[pops$genotipo %in% unique(df$genotipo), ]
pops <- unique(pops)
cat(sprintf("   %d individuos con población asignada\n", nrow(pops)))

# ── join ──────────────────────────────────────────────────────────────────────
df_long <- df_long |>
  left_join(pops, by = "id")

# ── codificar alelos (0/1/2 → diploid; -9 → missing) ─────────────────────────
df_long <- df_long |>
  mutate(
    allele1 = case_when(
      coding == "0"  ~ 1L,
      coding == "1"  ~ 1L,
      coding == "2"  ~ 2L,
      coding == "-9" ~ 0L
    ),
    allele2 = case_when(
      coding == "0"  ~ 1L,
      coding == "1"  ~ 2L,
      coding == "2"  ~ 2L,
      coding == "-9" ~ 0L
    )
  )

# ── filtrar loci con solo missing ─────────────────────────────────────────────
df_long <- df_long |>
  group_by(locID) |>
  filter(any(!(allele1 == 0L | allele2 == 0L))) |>
  ungroup()

# ── formato BA3 ───────────────────────────────────────────────────────────────
ba3 <- df_long |>
  select(
    indivID = genotipo,
    popID   = popname_str,
    locID,
    allele1,
    allele2
  ) |>
  mutate(popID = gsub("\\s+", "_", popID))

cat(sprintf("── Individuos: %d | Loci: %d | Poblaciones: %d\n",
            n_distinct(ba3$indivID),
            n_distinct(ba3$locID),
            n_distinct(ba3$popID)))
cat("   Poblaciones:", paste(unique(ba3$popID), collapse = ", "), "\n")

# ── exportar ──────────────────────────────────────────────────────────────────
write.table(ba3, file = outfile, sep = " ",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("✓ Guardado:", outfile, "\n")
