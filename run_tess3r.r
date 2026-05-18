#!/usr/bin/env Rapp
#| description: Corre tess3r sobre un archivo STR y guarda el objeto RDS resultante.

#| description: Ruta al archivo STR (tab-separated, con header)
str_file <- NULL

#| description: Ruta al CSV de coordenadas (con columnas id, lat, lon)
coord_file <- NULL

#| description: Rango de K a explorar, ej: '1:15'
k_range <- "1:15"

#| description: Método de optimización (qp o projected.ls)
method <- "qp"

#| description: Número de cores OpenMP
cores <- 6L

#| description: Nombre del archivo RDS de salida
output <- "tess3.rds"

library(tess3r)

# --- cargar datos ---
str_df <- read.csv(str_file, sep = "\t", header = TRUE)
str_df  <- str_df[seq(2, nrow(str_df), 2), ]

coord <- read.csv(coord_file, row.names = 1)
coord$pop <- NULL
coord <- unique(coord[coord$id %in% str_df$X, ])

namev            <- str_df$X
rownames(str_df) <- namev
str_df$X         <- NULL

rownames(coord) <- coord$id
coord$id        <- NULL
coord           <- coord[rownames(str_df), ]

# --- convertir formato ---
tss <- tess2tess3(
  dataframe    = str_df,
  TESS         = FALSE,
  diploid      = FALSE,
  FORMAT       = 1,
  extra.row    = 0,
  extra.column = 0
)

# --- correr tess3 ---
k_seq <- eval(parse(text = k_range))

tess3.obj <- tess3(
  X               = tss$X,
  coord           = tss$coord,
  K               = k_seq,
  method          = method,
  ploidy          = 1,
  openMP.core.num = cores
)

saveRDS(list(tess3 = tess3.obj, coord = coord), output)
cat("Objeto guardado en:", output, "\n")
