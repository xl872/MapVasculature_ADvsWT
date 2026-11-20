"""Reorient WT10 75um raw volume for SR: flips only, no resampling, orthogonal 0.075 mm affine.
Original header (sform_code 1): axcodes L,I,P (oblique). The other four raw volumes are R,I,A.
Action: flip array axis 0 (L->R) and axis 2 (P->A). No permutation (axis 1 is I in both).
New affine: orthogonal RIA, 0.075 mm, with the world position of the array centre kept
(centre voxel index (n-1)/2 is unchanged by flips). Data values and dtype unchanged."""
import sys, numpy as np, nibabel as nib
src, dst = sys.argv[1:3]
x = nib.load(src); d = x.dataobj.get_unscaled()   # stored int16 values, no scaling (scl_slope is NaN)
A0 = x.affine.copy()
d2 = np.ascontiguousarray(d[::-1, :, ::-1])
M = np.array([[0.075, 0, 0], [0, 0, 0.075], [0, -0.075, 0]])
cv = (np.array(d.shape) - 1) / 2.0
cw = A0[:3, :3] @ cv + A0[:3, 3]
A = np.eye(4); A[:3, :3] = M; A[:3, 3] = cw - M @ cv
h = x.header.copy(); h.set_data_dtype(d.dtype); h["scl_slope"] = np.nan; h["scl_inter"] = np.nan
img = nib.Nifti1Image(d2, A, h); img.set_sform(A, code=1); img.set_qform(A, code=0)
nib.save(img, dst)
np.set_printoptions(precision=6, suppress=True)
print('source', src, 'shape', d.shape, 'dtype', d.dtype, 'axcodes', nib.aff2axcodes(A0)); print(A0)
print('operation: data = data[::-1, :, ::-1]  (flip axis0 and axis2, no permutation, no resampling)')
print('new affine (sform_code 1, qform_code 0), axcodes', nib.aff2axcodes(A)); print(A)
print('world position of array centre kept:', cw)
r = nib.load(dst); assert r.get_data_dtype() == d.dtype and np.array_equal(r.dataobj.get_unscaled(), d2)
print('check: values identical after flip; sum', int(d.astype(np.int64).sum()), int(r.dataobj.get_unscaled().astype(np.int64).sum()))
