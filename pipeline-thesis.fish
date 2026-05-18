#!/usr/bin/env fish

# ─────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────
set SPECIES   "maulinus"
set RUN_DATE  (date +%Y%m%d)
set PFX       "{$SPECIES}_{$RUN_DATE}"
set WORKDIR   (pwd)
set RAPP      "Rapp"

# ─────────────────────────────────────────────
# INPUTS
# ─────────────────────────────────────────────
set RAW_DART   "$WORKDIR/raw_dart.csv"
set COV_FILE   "$WORKDIR/covariables.csv"
set ENV_FILE   "$WORKDIR/envar_final-corto.csv"
set COORD_FILE "$WORKDIR/coords.csv"

# ─────────────────────────────────────────────
# STEP 1 — filtrado
# ─────────────────────────────────────────────
echo "[1] Filtrado DArT..."
$RAPP $WORKDIR/lectura_et_filtrado.r \
    --filename    $RAW_DART \
    --covfilename $COV_FILE \
    --outname     "$WORKDIR/step1_gene/$PFX" \
    --savepath    "$WORKDIR/step1_gene"
or exit 1

set VCF    "$WORKDIR/step1_gene/{$PFX}.vcf"
set STR2   "$WORKDIR/step1_gene/{$PFX}.str"
set STR1   "$WORKDIR/step1_gene/{$PFX}_ploid1.str"
set GL_RDS "$WORKDIR/step1_obj/{$PFX}.rds"
set FASTA  "$WORKDIR/step1_gene/{$PFX}.fasta"

# ─────────────────────────────────────────────
# STEP 2 — TESS3r
# ─────────────────────────────────────────────
echo "[2] TESS3r..."
$RAPP $WORKDIR/run_tess3r.r \
    --str_file   $STR1 \
    --coord_file $COORD_FILE \
    --cores      6 \
    --output     "$WORKDIR/step2_obj/{$PFX}_tess3.rds"
or exit 1

$RAPP $WORKDIR/plot_tess3r.r \
    --rds_file   "$WORKDIR/step2_obj/{$PFX}_tess3.rds" \
    --coord_file $COORD_FILE \
    --outdir     "$WORKDIR/step2_figs"
or exit 1

# ─────────────────────────────────────────────
# STEP 3 — LFMM2
# ─────────────────────────────────────────────
echo "[3] LFMM2..."
$RAPP $WORKDIR/run_lfmm2.r \
    --str_file  $STR1 \
    --env_file  $ENV_FILE \
    --output    "$WORKDIR/step3_obj/{$PFX}_lfmm.rds" \
    --outdir    "$WORKDIR/step3_textfiles"
or exit 1

set VARLEA "$WORKDIR/step3_textfiles/varimportance_lea.csv"
set SNPS_SELECTIVE "$WORKDIR/step3_textfiles/sitios_selectivos.txt"

# ─────────────────────────────────────────────
# STEP 4 — GF + turnover
# ─────────────────────────────────────────────
echo "[4] Random Forest + Turnover..."
$RAPP $WORKDIR/random_forest_model.R \
    --vcf_file  $VCF \
    --env_file  $ENV_FILE \
    --output    "$WORKDIR/step4_obj/{$PFX}_forest_model.rds"
or exit 1

$RAPP $WORKDIR/composite_turnover_full.R \
    --forest_rds    "$WORKDIR/step4_obj/{$PFX}_forest_model.rds" \
    --env_file      $ENV_FILE \
    --outdir_text   "$WORKDIR/step4_textfiles" \
    --outdir_figs   "$WORKDIR/step4_figs" \
    --outdir_obj    "$WORKDIR/step4_obj"
or exit 1

set SNPS_OK      "$WORKDIR/step4_textfiles/snps_diagnostico_ok.csv"
set SNPS_PUTATIVE "$WORKDIR/step4_textfiles/snps_putativos.csv"
set VARIMGF      "$WORKDIR/step4_textfiles/var_importance_summary.csv"

# ─────────────────────────────────────────────
# STEP 5 — logical ops + format
# ─────────────────────────────────────────────
echo "[5] Clasificación loci + formateo..."

# 5a: logical ops → loci_clasificados.rds
$RAPP $WORKDIR/selective_loci_lfmm-gf.r \
    --lfmm_selective $SNPS_SELECTIVE \
    --gf_selective   $SNPS_OK \
    --gf_putative    $SNPS_PUTATIVE \
    --all_loci_rds   $GL_RDS \
    --gl_rds         $GL_RDS \
    --outfile        "$WORKDIR/step5_obj/{$PFX}_loci_clasificados.rds" \
    --outdir         "$WORKDIR/step5_textfiles"
or exit 1

set LOCI_RDS "$WORKDIR/step5_obj/{$PFX}_loci_clasificados.rds"

# 5b: TODO — script format_str_fasta.r (str + fasta simultáneo por categoría)
# $RAPP $WORKDIR/format_str_fasta.r \
#     --str_file  $STR2 \
#     --loci_rds  $LOCI_RDS \
#     --outdir    "$WORKDIR/step5_gene"

# 5c: TODO — script format_ba3.r
# $RAPP $WORKDIR/format_ba3.r \
#     --str_file  $STR2 \
#     --loci_rds  $LOCI_RDS \
#     --outdir    "$WORKDIR/step5_gene"

# ─────────────────────────────────────────────
# STEP 6 — IQ-TREE
# ─────────────────────────────────────────────
echo "[6] IQ-TREE por categoría..."
for CAT in true_neutral selected_and selected_or putative_local
    set FASTA_CAT "$WORKDIR/step5_gene/{$CAT}.fasta"
    set IQDIR     "$WORKDIR/step6_iqtree/iqtree-{$CAT}"
    mkdir -p $IQDIR
    iqtree2 -s $FASTA_CAT -T 4 \
        --prefix "$IQDIR/{$PFX}_{$CAT}" \
        -m GTR+G --quiet
    or exit 1
end

# ─────────────────────────────────────────────
# STEP 7 — multivariate
# ─────────────────────────────────────────────
echo "[7] Tests multivariados..."
$RAPP $WORKDIR/multivariate_comparation_tests.r \
    --genetic_distance_file "$WORKDIR/step6_iqtree/iqtree-true_neutral/{$PFX}_true_neutral.mldist" \
    --varimportance_lea     $VARLEA \
    --varimportance_gf      $VARIMGF \
    --env_file              $ENV_FILE \
    --coords_file           $COORD_FILE \
    --outfile               "$WORKDIR/step7_obj/{$PFX}_multivariate.rds" \
    --euler_png             "$WORKDIR/step7_figs/{$PFX}_euler.png"
or exit 1

# ─────────────────────────────────────────────
# STEP 8 — BayesAss
# ─────────────────────────────────────────────
echo "[8] BayesAss3 + plots..."
for CAT in true_neutral selected_and selected_or putative_local
    set STR_CAT "$WORKDIR/step5_gene/strfile_{$CAT}.str"
    BA3 -F $STR_CAT \
        -o "$WORKDIR/step8_textfiles/{$PFX}_{$CAT}_ba3out.txt"
    or exit 1
end

$RAPP $WORKDIR/plot_bayesass_panel.r \
    --main_title    "{$SPECIES} flujo génico" \
    --slots         4 \
    --title1        "Neutral" \
    --data1         "$WORKDIR/step8_textfiles/{$PFX}_true_neutral_ba3out.txt" \
    --title2        "Selectivo AND" \
    --data2         "$WORKDIR/step8_textfiles/{$PFX}_selected_and_ba3out.txt" \
    --title3        "Selectivo OR" \
    --data3         "$WORKDIR/step8_textfiles/{$PFX}_selected_or_ba3out.txt" \
    --title4        "Local" \
    --data4         "$WORKDIR/step8_textfiles/{$PFX}_putative_local_ba3out.txt" \
    --poptable_file $COORD_FILE \
    --outdir        "$WORKDIR/step8_figs"
or exit 1

echo "[DONE] Pipeline completo: $PFX"
