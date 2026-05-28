#!/usr/bin/env Rapp
#| description: Pipeline de filtrado SNP para datos DArT

#| description: Path al archivo SNP de DArT
filename <- ""

#| description: Path al archivo de covariables/poblaciones
covfilename <- ""

#| description: Keyword base para los archivos de salida
outname <- ""

#| description: folder para guardar los outputs
savepath <- ""

#| description: plink dir
plink<-""

# --- LIBS ---

library(dartRverse)
library(dartR.popgenomics)
library(dartR.base)
library(adegenet)
library(furrr)
library(parallel)
library(tidyverse)
library(stringr)


keyword <- outname
plink_path<-paste0(".",plink)
cat(sprintf("📂 filename:    %s\n", filename))
cat(sprintf("📂 covfilename: %s\n", covfilename))
cat(sprintf("🏷️  keyword:     %s\n", keyword))

# --- LECTURA ---
cat("\n📂 Leyendo datos DArT...\n")
gl <- gl.read.dart(filename = filename, covfilename = covfilename)
cat(sprintf("  gl: %d ind | %d loci | %d pops\n", nInd(gl), nLoc(gl), nPop(gl)))

# --- FILTROS SECUENCIALES ---
cat("\n🔍 Filtrando por read depth (5–75)...\n")
gl0 <- gl.filter.rdepth(gl, lower = 5, upper = 75)
cat(sprintf("  gl0: %d ind | %d loci (perdidos: %d)\n", nInd(gl0), nLoc(gl0), nLoc(gl)-nLoc(gl0)))

cat("\n🔍 Filtrando por MAF...\n")
gl1 <- gl.filter.maf(gl0, threshold = (1 / nInd(gl0)**2))
cat(sprintf("  gl1: %d ind | %d loci (perdidos: %d)\n", nInd(gl1), nLoc(gl1), nLoc(gl0)-nLoc(gl1)))

cat("\n🔍 Filtrando callrate por locus (threshold=0.5)...\n")
gl2 <- gl.filter.callrate(gl1, method = "loc", threshold = 0.5)
cat(sprintf("  gl2: %d ind | %d loci (perdidos: %d)\n", nInd(gl2), nLoc(gl2), nLoc(gl1)-nLoc(gl2)))

cat("\n🔍 Filtrando callrate por individuo (threshold=0.3)...\n")
gl3 <- gl.filter.callrate(gl2, method = "ind", threshold = 0.3)
cat(sprintf("  gl3: %d ind (perdidos: %d) | %d loci\n", nInd(gl3), nInd(gl2)-nInd(gl3), nLoc(gl3)))

cat("\n🔍 Filtrando secundarios...\n")
gl4 <- gl.filter.secondaries(gl3)
cat(sprintf("  gl4: %d ind | %d loci (perdidos: %d)\n", nInd(gl4), nLoc(gl4), nLoc(gl3)-nLoc(gl4)))

cat("\n🔍 Filtrando reproducibilidad (t=0.95)...\n")
gl5 <- gl.filter.reproducibility(gl4, t = 0.95)
cat(sprintf("  gl5: %d ind | %d loci (perdidos: %d)\n", nInd(gl5), nLoc(gl5), nLoc(gl4)-nLoc(gl5)))

cat("\n🔍 Filtrando monomorfos...\n")
gl6 <- gl.filter.monomorphs(gl5)
cat(sprintf("  gl6: %d ind | %d loci (perdidos: %d)\n", nInd(gl6), nLoc(gl6), nLoc(gl5)-nLoc(gl6)))

# --- CALLRATE DENTRO DE POPS ---
cat("\n🧬 Filtrando loci por callrate dentro de cada población...\n")
pops <- unique(as.vector(pop(gl6)))
cat(sprintf("  Poblaciones: %d → %s\n", length(pops), paste(pops, collapse=", ")))

loclists <- lapply(pops, \(p) {
  sub  <- gl6[pop(gl6)==p,]
  locs <- locNames(gl.filter.callrate(sub, threshold=0.5))
  cat(sprintf("    Pop '%s': %d loci pasan callrate\n", p, length(locs)))
  locs
})
counts <- sort(table(unlist(loclists)), decreasing=TRUE)
cat(sprintf("  Total loci únicos: %d | max count: %d | umbral: >= %d\n",
            length(counts), max(counts), max(counts)-2))

gl7 <- gl6[, locNames(gl6) %in% names(counts[counts >= (max(counts)-2)])]
cat(sprintf("  gl7: %d ind | %d loci\n", nInd(gl7), nLoc(gl7)))

# --- CALLRATE DURO ANTES DE LD ---
cat("\n🔍 Filtrando callrate global antes de LD (threshold=0.75)...\n")
gl7b <- gl.filter.callrate(gl7, method = "loc", threshold = 0.75)
cat(sprintf("  gl7b: %d ind | %d loci (perdidos: %d)\n", nInd(gl7b), nLoc(gl7b), nLoc(gl7)-nLoc(gl7b)))

# --- LD EN CHUNKS ---
chunk_size <- 175
step       <- 10
starts     <- seq(1, nLoc(gl7b) - chunk_size + 1, by = step)
ends       <- starts + chunk_size - 1
nbins      <- length(starts)

cat(sprintf("\n⚙️  LD chunking: %d loci | chunk=%d | step=%d | chunks=%d\n",
            nLoc(gl7b), chunk_size, step, nbins))

plan(multisession, workers = parallel::detectCores() - 1)
cat(sprintf("  Workers: %d\n", parallel::detectCores() - 1))

process_chunk <- function(i, gl_obj, starts, ends, nbins){
  tryCatch({
    library(adegenet)
    library(dartR.base)
    chunk_indices <- if (i == nbins) starts[i]:nLoc(gl_obj) else starts[i]:ends[i]
    values <- gl_obj[, chunk_indices]
    report <- gl.report.ld.map(values, ind.limit = 3, bins = 8)
    filt   <- gl.filter.ld(values, ld.report = report)
    filt$loc.names
  }, error = function(e){
    message(sprintf("💥 chunk fallido %d: %s", i, conditionMessage(e)))
    NULL
  })
}

cat("\n🚀 Ejecutando LD en paralelo...\n")
results <- future_map(
  seq_len(nbins),
  ~process_chunk(.x, gl7b, starts, ends, nbins),
  .options = furrr_options(
    globals  = list(gl_obj = gl7b, starts = starts, ends = ends, nbins = nbins, process_chunk = process_chunk),
    packages = c("adegenet", "dartR.base", "dartRverse"),
    seed     = TRUE
  )
)
cat(sprintf("  OK: %d | fallidos: %d\n", sum(!sapply(results, is.null)), sum(sapply(results, is.null))))

# --- CONSOLIDACIÓN ---
cat("\n📊 Consolidando resultados LD...\n")
vec_unico <- unlist(results)
cat(sprintf("  Entradas totales: %d\n", length(vec_unico)))

resultado <- tibble(valor = vec_unico) |>
  count(valor, name = "n") |>
  mutate(numberid = as.numeric(str_split_i(valor, "-", i = 1))) |>
  arrange(numberid) |>
  filter(n > 14)
cat(sprintf("  Loci retenidos (n>14 chunks): %d\n", nrow(resultado)))

gl8 <- gl7b[, resultado$valor]
cat(sprintf("  gl8: %d ind | %d loci\n", nInd(gl8), nLoc(gl8)))

# --- IMPUTACIÓN ---
cat("\n🔧 Imputando por frecuencia...\n")
gl10 <- gl.impute(x = gl8, method = "frequency")
cat(sprintf("  gl10 (final): %d ind | %d loci\n", nInd(gl10), nLoc(gl10)))

# --- EXPORTACIÓN ---
cat(sprintf("\n💾 Exportando con keyword '%s'...\n", keyword))

# --- EXPORTACIÓN ---
cat(sprintf("\n💾 Exportando con keyword '%s'...\n", keyword))
saveRDS(gl10, file = paste0(keyword, ".rds"))
gl2structure(gl10, outfile = paste0(basename(keyword), ".str"),        outpath = savepath, ploidy = 2)
gl2structure(gl10, outfile = paste0(basename(keyword), "_ploid1.str"), outpath = savepath, ploidy = 1)
gl2vcf(gl10,       outfile = basename(keyword),                        outpath = savepath, plink.bin.path = plink_path)
gl2fasta(gl10,     outfile = paste0(basename(keyword), ".fasta"),      outpath = savepath, method = 3)
gl2plink(gl10,     outfile = basename(keyword),                        outpath = savepath, plink.bin.path = plink_path, bed.files = TRUE)
gl2gds(gl10,       outfile = paste0(basename(keyword), ".gds"),        outpath = savepath)
cat("✅ Pipeline completo.\n")

cat("✅ Pipeline completo.\n")
