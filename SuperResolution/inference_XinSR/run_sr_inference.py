#!/usr/bin/env python
import os
import tensorflow as tf
import XinSR.utils
import XinSR.networks
import nibabel as nib
import numpy as np
from scipy.ndimage import zoom

# Configuration
# Edit DATA_DIR / SUBJECTS to point at the raw 75 um NIfTI files (Harvard Dataverse: in_vivo_75um/WT|AD)
DATA_DIR = '../../../HarvardDataverse/in_vivo_75um/WT'
SUBJECTS = ['01']   # subject IDs, e.g. ['01','04']; WT02/WT03 files have no _75um suffix
for mousenum in SUBJECTS:
    input_image = os.path.join(DATA_DIR, 'AZ_Vasculature_WT'+mousenum+'_75um.nii')
    print(input_image)
    load_weights = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'XinSR.h5')
    output_image = input_image[:-4] + '.sr.nii.gz'
    output_image_zoom = input_image[:-4] + '.linear.nii.gz'
    gpu = None

    # CNN architecture
    enc_nf = [16, 32, 32]
    dec_nf = [32, 32, 32, 32, 16, 16]
    factor = 2

    # set up devices
    tf.compat.v1.disable_eager_execution()
    device, nb_devices = XinSR.utils.setup_device(gpu)

    # prepare model and load weights
    model = XinSR.networks.Unet(inshape=[None, None, None, 1],
                                nb_features=[enc_nf, dec_nf],
                                name='unet')
    model.load_weights(load_weights)

    # Load image and blow up to make prediction
    print('Loading image')
    x = nib.load(input_image)
    volume = np.zeros([192,128,192])
    volume[0:x.shape[0], 0:x.shape[1], 0:x.shape[2]]=np.squeeze(x.get_fdata())
    volume[volume<0] = 0
    sorted = np.sort(volume, axis=None)
    idx = int(len(sorted) * 0.98)
    maxi = sorted[idx]
    volume = volume / maxi
    aff = x.affine
    print('Zooming in')
    zoomed = zoom(volume, factor, order=1, prefilter=False)

    # Prediction
    print('Predicting')
    pred = np.squeeze(model.predict(zoomed[np.newaxis,...,np.newaxis]))
    pred = maxi * (zoomed + pred)


    print('Writing to disk')
    header = nib.Nifti1Header()
    aff[0:-1, 0:-1] = aff[0:-1, 0:-1] / factor
    aff[0:-1, -1] = aff[0:-1, -1] - np.matmul(aff[0:-1, 0:-1], np.array([factor/2.0-0.5, factor/2.0-0.5, factor/2.0-0.5]))
    nifti = nib.Nifti1Image(pred, aff, header)
    nib.save(nifti, output_image)
    nifti = nib.Nifti1Image(zoomed * maxi, aff, header)
    nib.save(nifti, output_image_zoom)

    print('All done! freeview ' + input_image + ' ' + output_image_zoom + ' ' + output_image)

