#!/usr/bin/env python
"""SR of one volume with the computation of SuperResolution/inference_XinSR/run_sr_inference.py.

Same network config, padding to 192x128x192, 98th-percentile scaling, order-1 zoom
without prefilter, residual add and affine update as run_sr_inference.py; the input
and output paths come from the command line.

usage: run_sr_one.py <input.nii> <output.sr.nii.gz> <output.linear.nii.gz>
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SR_DIR = os.path.abspath(os.path.join(HERE, '..', 'SuperResolution', 'inference_XinSR'))
sys.path.insert(0, SR_DIR)

import numpy as np
import nibabel as nib
import tensorflow as tf
from scipy.ndimage import zoom
import XinSR.utils
import XinSR.networks

input_image, output_image, output_image_zoom = sys.argv[1:4]
load_weights = os.path.join(SR_DIR, 'XinSR.h5')
gpu = None

enc_nf = [16, 32, 32]
dec_nf = [32, 32, 32, 32, 16, 16]
factor = 2

tf.compat.v1.disable_eager_execution()
device, nb_devices = XinSR.utils.setup_device(gpu)
model = XinSR.networks.Unet(inshape=[None, None, None, 1],
                            nb_features=[enc_nf, dec_nf], name='unet')
model.load_weights(load_weights)

print('Loading image', input_image, flush=True)
x = nib.load(input_image)
volume = np.zeros([192, 128, 192])
volume[0:x.shape[0], 0:x.shape[1], 0:x.shape[2]] = np.squeeze(x.get_fdata())
volume[volume < 0] = 0
srt = np.sort(volume, axis=None)
idx = int(len(srt) * 0.98)
maxi = srt[idx]
volume = volume / maxi
aff = x.affine
print('input shape', x.shape, 'maxi', maxi, flush=True)
zoomed = zoom(volume, factor, order=1, prefilter=False)

print('Predicting', flush=True)
pred = np.squeeze(model.predict(zoomed[np.newaxis, ..., np.newaxis]))
pred = maxi * (zoomed + pred)

header = nib.Nifti1Header()
aff[0:-1, 0:-1] = aff[0:-1, 0:-1] / factor
aff[0:-1, -1] = aff[0:-1, -1] - np.matmul(aff[0:-1, 0:-1], np.array([factor/2.0-0.5]*3))
nib.save(nib.Nifti1Image(pred, aff, header), output_image)
nib.save(nib.Nifti1Image(zoomed * maxi, aff, header), output_image_zoom)
print('done', output_image, flush=True)
