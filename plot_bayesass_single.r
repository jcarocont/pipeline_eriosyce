#!/usr/bin/env Rapp
#| description: Genera un mapa de flujo migratorio desde un output BayesAss3

#| description: Título del gráfico y nombre base del PNG de salida
main_title <- NULL

#| description: Path al BA3out
data <- NULL

#| description: Path al CSV de tabla de poblaciones (con columnas: popname_str, lat, lon)
poptable_file <- NULL

#| description: Umbral mínimo de migración para dibujar flechas
mig_threshold <- 0.01

#| description: Directorio de salida
outdir <- "."

library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── validación ────────────────────────────────────────────────────────────────
if (is.null(main_title))   stop("--main_title es requerido")
if (is.null(data))         stop("--data es requerido")
if (is.null(poptable_file)) stop("--poptable_file es requerido")

# ── coordenadas ───────────────────────────────────────────────────────────────
coords_pop <- read.csv(poptable_file) |>
  group_by(popname_str) |>
  summarise(
    lat   = mean(lat, na.rm = TRUE),
    lon   = mean(lon, na.rm = TRUE),
    n_ind = n(),
    .groups = "drop"
  )

# ── parsear BA3out ────────────────────────────────────────────────────────────
out <- readLines(data, warn = FALSE)

idx <- grep("->", out)
lookup_df <- out[idx[2]] |>
  str_split(" ", simplify = TRUE) |>
  as.character() |>
  keep(~ nchar(trimws(.)) > 0) |>
  str_split("->", simplify = TRUE) |>
  as.data.frame(stringsAsFactors = FALSE) |>
  setNames(c("id", "population")) |>
  mutate(id = as.integer(id), population = gsub("_", " ", population)) |>
  as_tibble()

pop_lookup <- setNames(lookup_df$population, lookup_df$id)

start         <- grep("Migration",  out)
end           <- grep("Inbreeding", out)
migration_vec <- out[(start + 2):(end - 2)]

raw_matches <- migration_vec |>
  map(~ str_extract_all(
    .x, "m\\[(\\d+)\\]\\[(\\d+)\\]:\\s*([0-9.]+)\\(([0-9.]+)\\)"
  )[[1]]) |>
  unlist()

if (length(raw_matches) == 0)
  stop(sprintf("No se encontraron entradas m[i][j] en: %s", data))

df_long <- tibble(raw = raw_matches) |>
  mutate(
    row     = as.integer(str_match(raw, "m\\[(\\d+)\\]")[, 2]),
    col     = as.integer(str_match(raw, "m\\[\\d+\\]\\[(\\d+)\\]")[, 2]),
    value   = as.numeric(str_match(raw, ":\\s*([0-9.]+)")[, 2]),
    se      = as.numeric(str_match(raw, "\\(([0-9.]+)\\)")[, 2]),
    row_pop = pop_lookup[as.character(row)],
    col_pop = pop_lookup[as.character(col)]
  ) |>
  filter(!is.na(row), !is.na(col))

all_pops    <- lookup_df$population
df_complete <- expand_grid(row_pop = all_pops, col_pop = all_pops) |>
  left_join(df_long, by = c("row_pop", "col_pop"))

mig_matrix <- df_complete |>
  select(row_pop, col_pop, value) |>
  pivot_wider(names_from = col_pop, values_from = value, values_fill = NA) |>
  column_to_rownames("row_pop") |>
  as.matrix()

mig_df <- as.data.frame(mig_matrix) |>
  rownames_to_column("from") |>
  pivot_longer(-from, names_to = "to", values_to = "m") |>
  left_join(coords_pop, by = c("from" = "popname_str")) |>
  rename(lat_from = lat, lon_from = lon) |>
  left_join(coords_pop, by = c("to" = "popname_str")) |>
  rename(lat_to = lat, lon_to = lon)

self_flow <- mig_df |>
  filter(from == to) |>
  left_join(coords_pop, by = c("from" = "popname_str"))

mig_df <- mig_df |>
  filter(from != to, m >= mig_threshold, !is.na(lon_from), !is.na(lon_to))

if (nrow(mig_df) == 0)
  warning(sprintf("Sin flechas con mig_threshold = %g", mig_threshold))

# ── mapa base ─────────────────────────────────────────────────────────────────
chile <- ne_countries(scale = "medium", country = "Chile", returnclass = "sf")
pad   <- 1
bbox  <- st_bbox(c(
  xmin = min(coords_pop$lon) - pad, xmax = max(coords_pop$lon) + pad,
  ymin = min(coords_pop$lat) - pad, ymax = max(coords_pop$lat) + pad
), crs = st_crs(4326))
chile_crop <- st_crop(chile, bbox)

# ── plot ──────────────────────────────────────────────────────────────────────
m_range      <- range(mig_df$m,    na.rm = TRUE)
self_m_range <- range(self_flow$m, na.rm = TRUE)

p <- ggplot() +
  geom_sf(data = chile_crop, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_curve(
    data = mig_df |>
      mutate(
        m_scaled = (m - m_range[1]) / diff(m_range),
        m_vis    = pmin(pmax(m_scaled ^ 1.5, 0), 1)
      ) |>
      arrange(m),
    aes(lon_from, lat_from,
        xend = lon_to, yend = lat_to,
        colour = m, linewidth = m_vis, alpha = m_vis),
    curvature = 0.2,
    arrow = arrow(length = unit(0.01, "cm"))
  ) +
  geom_point(
    data = self_flow |> filter(!is.na(lon)),
    aes(lon, lat, size = m),
    shape = 21, fill = "red", color = "black", alpha = 0.7
  ) +
  scale_colour_viridis_c(
    limits = m_range, breaks = m_range,
    labels = scales::label_number(accuracy = 0.01),
    option = "turbo", name = "Migración"
  ) +
  scale_linewidth(range = c(0.05, 2), guide = "none") +
  scale_alpha(range = c(0, 1),        guide = "none") +
  scale_size(range = c(2, 8), limits = self_m_range, name = "Self-flow") +
  coord_sf(expand = FALSE) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  labs(title = main_title) +
  theme_minimal() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "left",
    axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1)
  )

# ── guardar ───────────────────────────────────────────────────────────────────
slug    <- gsub("[^a-zA-Z0-9]", "_", tolower(main_title))
outfile <- file.path(outdir, paste0(slug, ".png"))
ggsave(outfile, p, width = 7, height = 10, dpi = 300)
cat(sprintf("✓ Guardado: %s\n", outfile))

write.csv(mig_df,    file.path(outdir, paste0(slug, "_mig_df.csv")),    row.names = FALSE)
write.csv(self_flow, file.path(outdir, paste0(slug, "_self_flow.csv")), row.names = FALSE)
cat("✓ CSVs exportados\n")
