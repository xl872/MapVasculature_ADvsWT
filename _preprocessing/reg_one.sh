#!/bin/bash
# Step 3: register one SR volume to the WT06 SR reg1 base, then its raw 75um volume to that SR reg1.
# Options as in the AFNI history of the SR animals registered to the WT06 base (e.g. AD09, WT07):
#   3dUnifize <sr> -clfrac 0.5 ; 3dAllineate -base WT06_75um.sr_reg1+orig -source Unifized -autoweight -mast_dxyz 0.03 -lpa -twopass -cmass -final wsinc5
# regraw as in the history of e.g. AD09_75um.nii_regraw:
#   3dUnifize -clfrac 0.25 -prefix raw <raw75> ; 3dAllineate -base <own SR reg1> -source raw+orig -autoweight -mast_dxyz 0.03 -lpa -twopass -cmass -final wsinc5
# No 3drefit -orient: the SR inputs read as LSP in AFNI, like AD09 (registered without refit).
set -euo pipefail
ID=$1; SR=$2; RAW=$3
VESS=/path/to/vess   # <-- set this: WORK_DIR with AD/, WT/ (WT06 SR reg1 base in WT/)
BASE=$VESS/WT/AZ_Vasculature_WT06_75um.sr_reg1+orig
export AFNI_COMPRESSOR=NONE OMP_NUM_THREADS=16
D=$VESS/reg/$ID; mkdir -p $D; cd $D   # output folder per animal
cp $SR AZ_Vasculature_${ID}_75um.sr.nii.gz
cp $RAW AZ_Vasculature_${ID}_75um.nii
3dUnifize AZ_Vasculature_${ID}_75um.sr.nii.gz -clfrac 0.5
SRC=$(ls Unifized+*.HEAD | head -1); SRC=${SRC%.HEAD}
3dAllineate -base $BASE -source $SRC -autoweight -mast_dxyz 0.03 -lpa -twopass -cmass -final wsinc5 -prefix AZ_Vasculature_${ID}_75um.sr.nii.gz_reg1 -1Dmatrix_save ${ID}_sr_reg1.aff12.1D
3dUnifize -clfrac 0.25 -prefix raw AZ_Vasculature_${ID}_75um.nii
3dAllineate -base AZ_Vasculature_${ID}_75um.sr.nii.gz_reg1+orig -source raw+orig -autoweight -mast_dxyz 0.03 -lpa -twopass -cmass -final wsinc5 -prefix AZ_Vasculature_${ID}_75um.nii_regraw -1Dmatrix_save ${ID}_regraw.aff12.1D
echo DONE $ID
