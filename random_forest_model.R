#!/usr/bin/env Rapp
#| description: Path al VCF de entrada
vcf_file <- NULL

#| description: Path al CSV de variables ambientales
env_file <- NULL

#| description: Path de salida para el RDS de modelos
output <- "forest_model.rds"

library(ranger)
library(vcfR)
library(adegenet)
library(ape)
library(LEA)
library(tidyverse)
library(furrr)
library(rlang)

train_snp_model <- function(snp_vector, env, ntree=500, na_act="na.learn") {
  
  data_train <- cbind(snp=snp_vector, env)
  
  formula <- as.formula(paste0("snp", " ~ ."))
  
  tryCatch({
    ranger(
      formula,
      data = data_train,
      num.trees = ntree,
      mtry = 2,
      keep.forest = TRUE,
      write.forest = TRUE,
      min.node.size = 2,
      importance = "impurity",
      keep.inbag = TRUE,
      na.action = na_act,
      save.memory = TRUE,
      seed=TRUE,
      num.threads = 1
    )
  }, error = function(e) {
    warning(paste("Error entrenando modelo:", e$message))
    return(e$message)
  })
}

plan(multicore, workers = 8)
options(future.globals.maxSize = 5e9,future.rng.onMisuse = "ignore",future.stdout = TRUE)

vcf <- read.vcfR(vcf_file)

genlight_obj <- vcfR2genlight(vcf)

ind_names <- indNames(genlight_obj)|>
  str_split_fixed("genotype", 2) |>
  (\(x) paste0("genotype", unique(x[,2])))()

message("creating struct-like df")

df_structure <- as.data.frame(genlight_obj)
rownames(df_structure)<-ind_names
message("freeing memory")
rm(vcf); gc()
rm(genlight_obj); gc()
message("loading env matix")

env <- read.csv(env_file) |>
  column_to_rownames("id") |>
  select(-2,-X) |>
  as.matrix()

env <- env[ind_names,]
message("running models")



chunk_size <- 1000
n_snps <- ncol(df_structure)
lista <- split(1:n_snps, ceiling(seq_along(1:n_snps)/chunk_size))

all_models <- list()
last_i <- lista[[length(lista)]]


for (i in lista) {
  mess_str<-paste0("computando modelos para snps entre ", min(i), " - ", max(i))
  message(mess_str)
  models <- future_map(colnames(df_structure)[i], \(s) {
    x <- df_structure[[s]][!is.na(df_structure[[s]])]
    filt <- rownames(df_structure)[!is.na(df_structure[[s]])]
    env_f <- env[filt, , drop = FALSE]
    model <- train_snp_model(x, env_f, na_act = "na.omit")
    message("modelo entrenado exitosamente, para el SNP : ", s)
    setNames(list(list(model = model, percent = length(x)/length(df_structure[[s]]))), s)
  })
  models<-unlist(models, recursive = FALSE)
  if (identical(i, last_i)) {
        message("¡ÚLTIMO BLOQUE PROCESADO! guardando el chunck")
        saveRDS(models, "forest_model_660chunck.rds")
    }
  all_models <- c(all_models, models)
  rm(models);gc()
}


message("models finished")


message("saving models")
saveRDS(all_models, "forest_model.rds")
message("modelos entrenados y guardados. Ya hicimos lo mas complicado, ahora viene lo mas dificil")
