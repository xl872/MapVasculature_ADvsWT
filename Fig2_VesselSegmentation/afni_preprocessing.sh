#!/usr/bin/env bash
# AFNI preprocessing used for the vessel analysis (Fig 2A, Methods 2.5).
# The commands and options are those stored in the AFNI header history
# (HISTORY_NOTE) of each registered dataset. AFNI version: AFNI_23.0.02
# (linux_centos_7_64).
#
# Inputs : AZ_Vasculature_<ID>_sr.nii.gz  (CNN-SR output, 37.5 um)
#          AZ_Vasculature_<ID>_75um.nii   (raw 75 um)
#          ambmc_30+orig                  (AMBMC 30 um template, see HarvardDataverse/)
# Outputs: <ID>_sr_reg1+orig              -> myseg_tissue_top2to50.m
#          <ID>_75um.nii_regraw+orig      -> raw intensity used by the group scripts
#
# Per-subject settings recorded in the headers:
#   * 3dUnifize -clfrac: default for WT01-06, AD02-08, AD07_1/2;
#                        0.5 for WT07, WT09, WT10, WT12, AD09, AD10, AD11; 0.25 for WT11.
#   * 3dAllineate -base: ambmc_30+orig for most subjects; the registered WT06 SR
#     volume WT06_75um.sr_reg1+orig for WT07, WT09, WT10, WT12, AD09, AD10, AD11.
#   * WT01 and AD08 (both "_sr_distort") were registered twice (regtempt -> reg1).
#   * 3drefit before unifize: AD02 '-orient RPI' then '-orient LSA';
#     WT05/WT07 '-xdel 0.038 -ydel 0.038 -zdel 0.038'; WT07 '-duporigin WT06'.
#   * AD02-AD07_2: 3drefit -orient LSA Unifized+tlrc after unifize.
set -euo pipefail

TEMPLATE=ambmc_30+orig

# ---- 1. SR volume: bias-field correction + registration to AMBMC 30 um ----
reg_sr() {          # $1 = SR input file, $2 = output prefix, $3 = base (default template)
  local in=$1 out=$2 base=${3:-$TEMPLATE}
  3dUnifize -prefix Unifized "$in"                      # (add -clfrac per subject, see above)
  3drefit -orient LSA Unifized+tlrc
  3dAllineate -base "$base" -source Unifized+tlrc \
      -mast_dxyz 0.03 -autoweight -lpa -twopass -cmass -final wsinc5 \
      -prefix "$out"
  rm -f Unifized+tlrc.*
}

# ---- 2. raw 75 um volume: unifize + register to the subject's own SR reg1 ----
reg_raw() {         # $1 = ID  (e.g. AD02)
  local d=AZ_Vasculature_$1_75um
  3dUnifize -clfrac 0.25 -prefix raw ${d}.nii
  3drefit -orient LSA raw+orig
  3dAllineate -base AZ_Vasculature_$1_sr_reg1+orig -source raw+orig \
      -autoweight -mast_dxyz 0.03 -lpa -twopass -cmass -final wsinc5 \
      -prefix ${d}.nii_regraw
  rm -f raw+orig.*
}

# ---- 3. Fig S2/S3 only (AD02-AD05): raw on its own
#         75 um grid (regonly) and nearest-neighbour onto the 30 um grid (regNN)
reg_raw_native() {  # $1 = ID
  local d=AZ_Vasculature_$1_75um
  3dUnifize -clfrac 0.25 -prefix raw ${d}.nii
  3drefit -orient LSA raw+orig
  3dAllineate -base AZ_Vasculature_$1_sr_reg1+orig -source raw+orig -autoweight \
      -mast_dxyz 0.075 -lpa -twopass -cmass -final wsinc5 \
      -prefix ${d}_regonly -1Dmatrix_save ${d}_regonly.aff12.1D
  3dAllineate -1Dmatrix_apply ${d}_regonly.aff12.1D -source raw+orig \
      -master AZ_Vasculature_AD02_sr_reg1+orig -final NN -prefix ${d}_regNN
  rm -f raw+orig.*
}

# Example (AD group):
# for i in 02 03 04 05 06; do reg_sr AZ_Vasculature_AD${i}_sr.nii.gz AZ_Vasculature_AD${i}_sr_reg1; reg_raw AD${i}; done
# for i in 02 03 04 05;    do reg_raw_native AD${i}; done
