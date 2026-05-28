#!/usr/bin/env Rapp
#| description: Genera gráficos ggplot a partir del RDS de tess3r.

#| description: Ruta al RDS generado por run_tess3r
rds_file <- NULL

#| description: Ruta al CSV de coordenadas (con columnas id, lat, lon)
coord_file <- NULL

#| description: K a usar para barplot y mapa. Si es 0, se elige por mínimo de cross-entropy.
k <- 0L

#| description: Directorio de salida para los PNGs
outdir <- "."

#| description: Resolución de la interpolación del mapa (nro de puntos por eje)
map_res <- 300L

#| description: Theta (rango espacial) para la interpolación Krig
krig_theta <- 10

library(tess3r)
library(tidyverse)
library(ggplot2)
library(viridis)
library(patchwork)
library(fields)
library(sf)
library(rnaturalearth)
library(sp)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── cargar datos ────────────────────────────────────────────────────────────────
rds <- readRDS(rds_file)

# soporta tanto el formato nuevo (lista con $tess3 y $coord)
# como el formato viejo (objeto tess3 directo)
if (is.list(rds) && !is.null(rds$tess3)) {
  tess3.obj <- rds$tess3
  coord     <- rds$coord
} else {
  tess3.obj <- rds
  if (is.null(coord_file)) stop("Este RDS no incluye coordenadas. Usa --coord_file.")
  coord <- read.csv(coord_file, row.names = 1)
  if ("pop" %in% colnames(coord)) coord$pop <- NULL
  if ("id"  %in% colnames(coord)) { rownames(coord) <- coord$id; coord$id <- NULL }
}

# siempre se puede sobreescribir coord con --coord_file
if (!is.null(coord_file) && coord_file != "") {
  coord <- read.csv(coord_file, row.names = 1)
  if ("pop" %in% colnames(coord)) coord$pop <- NULL
  if ("id"  %in% colnames(coord)) { rownames(coord) <- coord$id; coord$id <- NULL }
}

# ── métricas por K ──────────────────────────────────────────────────────────────
output_list <- lapply(tess3.obj, function(i) {
  gif_val <- tryCatch(i$tess3.run[[1]]$gif, error = function(e) NA_real_)
  c(K = i$K, cross = i$crossentropy, rmse = i$rmse, gif = gif_val)
})
metrics <- do.call(rbind, output_list) |>
  as.data.frame() |>
  arrange(K)

# ── elegir K ────────────────────────────────────────────────────────────────────
k_opt <- if (k == 0L) {
  metrics$K[which.min(metrics$cross)]
} else {
  as.integer(k)
}
cat("K usado:", k_opt, "\n")

# ── paleta (hasta 15 clusters) ──────────────────────────────────────────────────
palette_base <- c(
  "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
  "#ffff33", "#a65628", "#f781bf", "#999999", "#66c2a5",
  "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f"
)

# ══════════════════════════════════════════════════════════════════════════════
# GRÁFICO 1 — Curva cross-entropy
# ══════════════════════════════════════════════════════════════════════════════
p_ce <- ggplot(metrics, aes(x = K, y = cross)) +
  geom_line(linewidth = 0.8, color = "grey40") +
  geom_point(size = 2.5, color = "grey20") +
  geom_point(
    data = filter(metrics, K == k_opt),
    aes(x = K, y = cross),
    color = "red", size = 4, shape = 18
  ) +
  annotate(
    "text",
    x = k_opt + 0.3, y = metrics$cross[metrics$K == k_opt],
    label = paste("K =", k_opt),
    color = "red", hjust = 0, size = 3.5
  ) +
  scale_x_continuous(breaks = metrics$K) +
  labs(
    x     = "Número de clusters (K)",
    y     = "Cross-entropy",
    title = "Selección de K óptimo"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(outdir, "cross_entropy.png"),
  p_ce, width = 7, height = 4.5, dpi = 300
)
cat("Guardado: cross_entropy.png\n")

# ══════════════════════════════════════════════════════════════════════════════
# GRÁFICO 2 — Barplot de ancestría
# ══════════════════════════════════════════════════════════════════════════════
q.matrix <- qmatrix(tess3.obj, K = k_opt)
colnames(q.matrix) <- paste0("K_", seq_len(ncol(q.matrix)))

admix_df <- cbind.data.frame(coord, q.matrix)
admix_df$id <- rownames(admix_df)
rownames(admix_df) <- NULL

cols <- palette_base[seq_len(k_opt)]

df_long <- admix_df |>
  pivot_longer(cols = starts_with("K_"), names_to = "cluster", values_to = "q") |>
  group_by(id) |>
  mutate(lat_id = first(lat)) |>
  ungroup() |>
  mutate(id = factor(id, levels = unique(id[order(-lat_id)])))

p_bar <- ggplot(df_long, aes(x = id, y = q, fill = cluster)) +
  geom_col(width = 1, color = "grey90", linewidth = 0.15) +
  scale_fill_manual(values = cols, name = "Cluster") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x     = "Individuos (N → S)",
    y     = "Proporción de ancestría",
    title = paste0("Ancestría por individuo (K = ", k_opt, ")")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
    axis.ticks.x    = element_blank(),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

ggsave(
  file.path(outdir, paste0("barplot_K", k_opt, ".png")),
  p_bar, width = 13, height = 5, dpi = 300
)
cat("Guardado: barplot_K", k_opt, ".png\n", sep = "")

# ══════════════════════════════════════════════════════════════════════════════
# GRÁFICO 3 — Mapa interpolado (Krig)
# ══════════════════════════════════════════════════════════════════════════════
coords_mat <- admix_df |> select(lon, lat) |> as.matrix()
Q_mat      <- admix_df |> select(starts_with("K_")) |> as.matrix()
K_n        <- ncol(Q_mat)

xg   <- seq(min(admix_df$lon), max(admix_df$lon), length.out = map_res)
yg   <- seq(min(admix_df$lat), max(admix_df$lat), length.out = map_res)
grid <- expand.grid(lon = xg, lat = yg)

interp_list <- lapply(seq_len(K_n), function(k_i) {
  Krig(x = coords_mat, Y = Q_mat[, k_i], theta = krig_theta)
})

Q_interp <- sapply(seq_len(K_n), function(k_i) predict(interp_list[[k_i]], grid))
colnames(Q_interp) <- paste0("K_", seq_len(K_n))

map_df <- grid |>
  bind_cols(as_tibble(Q_interp)) |>
  mutate(
    max_cluster = paste0("K_", max.col(across(starts_with("K_")), ties.method = "first"))
  )

# recortar al contorno de Chile
chile <- ne_countries(scale = "medium", country = "Chile", returnclass = "sf")

pts_sf  <- st_as_sf(admix_df, coords = c("lon", "lat"), crs = 4326)
pts_utm <- st_transform(pts_sf, 32719)

hull_ll <- st_transform(
  st_buffer(st_convex_hull(st_union(pts_utm)), 100000),
  4326
)
chile_cut <- st_intersection(chile, hull_ll)

chile_poly <- chile_cut |>
  st_union() |>
  st_cast("MULTIPOLYGON") |>
  st_coordinates()

inside <- point.in.polygon(
  point.x = map_df$lon,
  point.y = map_df$lat,
  pol.x   = chile_poly[, "X"],
  pol.y   = chile_poly[, "Y"]
) > 0

map_df_cut <- map_df[inside, ]

bbox_ll <- st_transform(
  st_buffer(st_as_sfc(st_bbox(pts_utm)), 10000),
  4326
)

p_map <- ggplot() +
  geom_sf(data = chile_cut, fill = "grey92", color = "grey50", linewidth = 0.3) +
  geom_raster(
    data  = map_df_cut,
    aes(x = lon, y = lat, fill = max_cluster),
    alpha = 0.85
  ) +
  geom_point(
    data  = admix_df,
    aes(x = lon, y = lat),
    color = "black", size = 0.8, shape = 16
  ) +
  scale_fill_manual(values = setNames(cols, paste0("K_", seq_len(k_opt))), name = "Cluster") +
  coord_sf(
    xlim = st_bbox(bbox_ll)[c("xmin", "xmax")],
    ylim = st_bbox(bbox_ll)[c("ymin", "ymax")]
  ) +
  labs(
    x     = "Longitud",
    y     = "Latitud",
    title = paste0("Ancestría interpolada (K = ", k_opt, ")")
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_line(color = "grey80", linewidth = 0.2))

ggsave(
  file.path(outdir, paste0("mapa_K", k_opt, ".png")),
  p_map, width = 7, height = 14, dpi = 300
)
cat("Guardado: mapa_K", k_opt, ".png\n", sep = "")
