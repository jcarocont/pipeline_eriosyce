#!/usr/bin/env Rapp
#| description: Genera un panel comparativo de mapas de flujo migratorio desde outputs BayesAss3

#| description: Título del panel y nombre base del PNG de salida
main_title <- NULL

#| description: Número de datasets a comparar (máximo 4)
slots <- 2L

#| description: Título subgráfico 1
title1 <- NULL
#| description: Path BA3out dataset 1
data1 <- NULL

#| description: Título subgráfico 2
title2 <- NULL
#| description: Path BA3out dataset 2
data2 <- NULL

#| description: Título subgráfico 3
title3 <- NULL
#| description: Path BA3out dataset 3
data3 <- NULL

#| description: Título subgráfico 4
title4 <- NULL
#| description: Path BA3out dataset 4
data4 <- NULL

#| description: Path al CSV de tabla de poblaciones (con columnas: popname_str, lat, lon)
poptable_file <- NULL

#| description: Umbral mínimo de migración para dibujar flechas
mig_threshold <- 0.01

#| description: Directorio de salida
outdir <- "."

library(tidyverse)
library(patchwork)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── validación ────────────────────────────────────────────────────────────────
if (is.null(main_title))   stop("--main_title es requerido")
if (is.null(poptable_file)) stop("--poptable_file es requerido")

all_titles <- list(title1, title2, title3, title4)
all_data   <- list(data1,  data2,  data3,  data4)

if (slots < 1 || slots > 4)
  stop("--slots debe estar entre 1 y 4")

missing_titles <- which(sapply(all_titles[1:slots], is.null))
missing_data   <- which(sapply(all_data[1:slots],   is.null))

if (length(missing_titles) > 0)
  stop(sprintf("--slots %d pero faltan títulos: %s",
               slots, paste0("title", missing_titles, collapse = ", ")))

if (length(missing_data) > 0)
  stop(sprintf("--slots %d pero faltan datos: %s",
               slots, paste0("data", missing_data, collapse = ", ")))

datasets <- map(1:slots, ~ list(
  title = all_titles[[.x]],
  path  = all_data[[.x]]
))

# ── coordenadas ───────────────────────────────────────────────────────────────
coords_pop <- read.csv(poptable_file) |>
  group_by(popname_str) |>
  summarise(
    lat   = mean(lat, na.rm = TRUE),
    lon   = mean(lon, na.rm = TRUE),
    n_ind = n(),
    .groups = "drop"
  )

# ── helper: parsear BA3out ────────────────────────────────────────────────────
parse_ba3 <- function(path, coords_pop) {
  out <- readLines(path, warn = FALSE)

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
    stop(sprintf("No se encontraron entradas m[i][j] en: %s", path))

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

  list(mig_df = mig_df, self_flow = self_flow, mig_matrix = mig_matrix)
}

# ── mapa base ─────────────────────────────────────────────────────────────────
chile <- ne_countries(scale = "medium", country = "Chile", returnclass = "sf")
pad   <- 1
bbox  <- st_bbox(c(
  xmin = min(coords_pop$lon) - pad, xmax = max(coords_pop$lon) + pad,
  ymin = min(coords_pop$lat) - pad, ymax = max(coords_pop$lat) + pad
), crs = st_crs(4326))
chile_crop <- st_crop(chile, bbox)

# ── función plot ──────────────────────────────────────────────────────────────
plot_mig_map <- function(mig_df, self_flow, titulo) {
  if (nrow(mig_df) == 0)
    warning(sprintf("Sin flechas para '%s' con mig_threshold = %g", titulo, mig_threshold))

  m_range      <- range(mig_df$m,    na.rm = TRUE)
  self_m_range <- range(self_flow$m, na.rm = TRUE)

  ggplot() +
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
    labs(title = titulo) +
    theme_minimal() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      legend.position = "left",
      axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1)
    )
}

# ── parsear y plotear ─────────────────────────────────────────────────────────
parsed <- map(datasets, ~ {
  cat(sprintf("Parseando: %s\n", .x$path))
  parse_ba3(.x$path, coords_pop)
})

plots <- map2(parsed, datasets, ~
  plot_mig_map(.x$mig_df, .x$self_flow, .y$title)
)

# ── panel con título global ───────────────────────────────────────────────────
panel <- wrap_plots(plots, nrow = 1, guides = "collect") +
  plot_annotation(
    title = main_title,
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  ) &
  theme(legend.position = "left")

slug    <- gsub("[^a-zA-Z0-9]", "_", tolower(main_title))
outfile <- file.path(outdir, paste0(slug, ".png"))
ggsave(outfile, panel, width = 5 * slots, height = 10, dpi = 300)
cat(sprintf("✓ Panel guardado: %s\n", outfile))

# ── exportar CSVs ─────────────────────────────────────────────────────────────
walk2(parsed, datasets, ~ {
  sub_slug <- gsub("[^a-zA-Z0-9]", "_", tolower(.y$title))
  write.csv(.x$mig_df,    file.path(outdir, paste0(slug, "_", sub_slug, "_mig_df.csv")),    row.names = FALSE)
  write.csv(.x$self_flow, file.path(outdir, paste0(slug, "_", sub_slug, "_self_flow.csv")), row.names = FALSE)
})
cat("✓ CSVs exportados\n")
