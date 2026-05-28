#!/usr/bin/env Rapp
#| description: Clasificación de loci selectivos vs neutrales (LFMM + GF) y exportación a FASTA via dartRverse.

#| description: Path al CSV de genes y cargas (LFMM)
lfmm_selective <- NULL

#| description: Path al CSV de sitios selectivos (GF)
gf_selective <- NULL

#| description: Path al CSV de sitios putativos (GF subset)
gf_subset <- NULL

#| description: Path al RDS del objeto genlight (dartRverse)
gl_rds <- NULL

#| description: Path de salida para el RDS de resultados
outfile_rds <- "loci_clasificados.rds"

#| description: Directorio de salida para los FASTA
outdir <- "."

#| description: Método gl2fasta (1=ambiguity codes, 2=random allele, 3=maj allele, 4=random haplotype)
fasta_method <- 1L

library(dartR.base)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# 1. CARGAR LOCI
# ══════════════════════════════════════════════════════════════════════════════
cat("── Cargando loci LFMM:", lfmm_selective, "\n")
genes     <- read.csv(lfmm_selective)
genes1    <- as.vector(genes)
nombres_genes <- names(genes1[genes1 > 0])
valnames  <- gsub("X", "", nombres_genes)
valnames  <- gsub("Response.", "", valnames)

selec_lmmf <- gsub(
  pattern     = "(\\d+)\\.(\\d+)\\.(\\w)\\.(\\w)",
  replacement = "\\1-\\2-\\3/\\4",
  x           = valnames
)
cat("   Loci LFMM crudos:", length(selec_lmmf), "\n")

cat("── Cargando loci GF:", gf_selective, "\n")
selec_rf <- as.vector(read.csv(gf_selective))[[1]]
cat("   Loci GF crudos:", length(selec_rf), "\n")
gl <- readRDS(gl_rds)
cat("── Cargando all_loci:", gl_rds, "\n")
all_loci <- locNames(gl)
cat("   Total loci de referencia:", length(all_loci), "\n")
gf_subset_names<-read.csv(gf_subset)
gf_subset_names<-unique(gf_subset_names$SNP_ID)
# ══════════════════════════════════════════════════════════════════════════════
# 2. INTERSECTAR CON all_loci Y CLASIFICAR
# ══════════════════════════════════════════════════════════════════════════════
cat("── Clasificando loci...\n")

selec_lmmf <- intersect(selec_lmmf, all_loci)
selec_rf   <- intersect(selec_rf,   all_loci)

cat("   Loci LFMM (tras intersección):", length(selec_lmmf), "\n")
cat("   Loci GF   (tras intersección):", length(selec_rf),   "\n")

selected_or  <- base::union(selec_rf, selec_lmmf)
selected_and <- intersect(selec_rf, selec_lmmf)
selected_xor <- unique(c(setdiff(selec_rf, selec_lmmf), setdiff(selec_lmmf, selec_rf)))
true_neutral <- setdiff(all_loci, selected_or)
neutral_rf   <- setdiff(all_loci, selec_rf)
neutral_lmmf <- setdiff(all_loci, selec_lmmf)

cat("   Selectivos solo GF   (xor):", length(setdiff(selec_rf, selec_lmmf)), "\n")
cat("   Selectivos solo LFMM (xor):", length(setdiff(selec_lmmf, selec_rf)), "\n")
cat("   Selectivos ambos (and)    :", length(selected_and), "\n")
cat("   Selectivos union (or)     :", length(selected_or),  "\n")
cat("   Neutrales verdaderos      :", length(true_neutral), "\n")

res <- list(
  selec_lmmf   = selec_lmmf,
  selec_rf     = selec_rf,
  selected_and = selected_and,
  selected_or  = selected_or,
  selected_xor = selected_xor,
  true_neutral = true_neutral,
  neutral_rf   = neutral_rf,
  neutral_lmmf = neutral_lmmf,
  putative_local=gf_subset_names
)

saveRDS(res, outfile)
cat("✓ RDS de clasificación guardado en:", outfile, "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# 3. EXPORTAR FASTA VIA DARTRVERSE
# ══════════════════════════════════════════════════════════════════════════════
cat("── Cargando genlight:", gl_rds, "\n")
cat("   Individuos:", nInd(gl), "| Loci:", nLoc(gl), "| Pops:", nPop(gl), "\n")
cat("   Poblaciones:", paste(popNames(gl), collapse = ", "), "\n")

# helper: filtrar → fasta con prints informativos
export_fasta <- function(gl, loci_vec, label) {
  cat("\n── Exportando FASTA:", label, "\n")

  # loci disponibles en el GL que coinciden con el subset
  loci_en_gl <- intersect(loci_vec, locNames(gl))
  cat("   Loci solicitados:", length(loci_vec),
      "| Encontrados en GL:", length(loci_en_gl), "\n")

  if (length(loci_en_gl) == 0) {
    cat("   ⚠ Sin loci comunes, se omite este FASTA.\n")
    return(invisible(NULL))
  }

  gl_sub <- gl.keep.loc(gl, loc.list = loci_en_gl, verbose = 0)
  cat("   GL filtrado → loci:", nLoc(gl_sub), "| ind:", nInd(gl_sub), "\n")

  fname <- paste0(label, ".fasta")
  gl2fasta(
    gl_sub,
    method  = fasta_method,
    outfile = fname,
    outpath = outdir,
    verbose = 0
  )
  cat("   ✓ Guardado:", file.path(outdir, fname), "\n")
}

export_fasta(gl, true_neutral, "true_neutral")
export_fasta(gl, gf_subset,"putative_local")
export_fasta(gl, selec_rf,    "selec_gf")
export_fasta(gl, selec_lmmf,  "selec_lfmm")

for (cat in names(res)) {
  gl_sub <- gl.keep.loc(gl, loc.list = intersect(res[[cat]], locNames(gl)), verbose = 0)
  gl2structure(
    gl_sub,
    outfile = paste0("strfile_", cat, ".str"),
    outpath = outdir,
    ploidy  = 2
  )
  cat(sprintf("✓ STR exportado: strfile_%s.str (%d loci)\n", cat, nLoc(gl_sub)))
}

cat("\n✓ Todo listo.\n")
