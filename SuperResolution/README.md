# Super-resolution (75 µm → 37.5 µm)

- `inference_XinSR/run_sr_inference.py`: makes `*.sr.nii.gz` (SR) and `*.linear.nii.gz` (trilinear) from a raw 75 µm NIfTI, using the trained weights `XinSR.h5`. **The SR volumes used in the paper are made with this script and these weights.** The network (`XinSR/networks.py`) is a 3D U-Net (TensorFlow/Keras; encoder 16-32-32, decoder 32-32-32, then 32-16-16 at full resolution; about 0.24 M parameters) that predicts a residual added to the trilinearly upsampled input. `XinSR.h5` records training with Adam (learning rate 1e-4), 700,000 iterations, 64³ patches. Set `DATA_DIR` and `SUBJECTS`. The input must be reoriented first so that its affine is diagonal.
- `training_ResSR/`: a separate PyTorch implementation of the training procedure with an EDSR-style residual network. It is **not** the code that produced `XinSR.h5` (the training script for the U-Net is not included). Its on-the-fly degradation of the ex vivo 37.5 µm volumes is in `ResSR/generators.py` and its network is in `ResSR/models.py`. Run `test_training.py`; `environment_training.yaml` is the conda environment.

All 22 animals are analysed from SR volumes made with `run_sr_inference.py` and `XinSR.h5`. `../_preprocessing/run_sr_one.py` is the same computation with the input and output paths from the command line; WT10 is reoriented first with `../_preprocessing/reorient_WT10.py`.

Training data: `ex_vivo_37.5um_training/` on Harvard Dataverse.
