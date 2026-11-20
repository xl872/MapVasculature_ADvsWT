#!/bin/bash

for i in 02 03 04 05
do
data=AZ_Vasculature_AD${i}_75um

3dUnifize -clfrac 0.25 -prefix raw ${data}.nii
3drefit -orient LSA raw+orig
3dAllineate -base AZ_Vasculature_AD${i}_sr_reg1+orig -source raw+orig -autoweight -mast_dxyz 0.075 -lpa -twopass -cmass -final wsinc5 -prefix ${data}_regonly -1Dmatrix_save ${data}_regonly.aff12.1D 

3dAllineate -1Dmatrix_apply ${data}_regonly.aff12.1D  \
    -source raw+orig \
    -master AZ_Vasculature_AD02_sr_reg1+orig \
    -final NN -prefix ${data}_regNN
rm raw+orig*

done
