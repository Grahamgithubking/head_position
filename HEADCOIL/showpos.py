import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from nibabel.affines import apply_affine

def show_slices(slices):
    """ Function to display row of image slices """
    fig, axes = plt.subplots(1, len(slices))
    for i, slice in enumerate(slices):
        axes[i].imshow(slice.T, cmap="gray", origin="lower", vmin=0, vmax=100)
        axes[i].set_axis_off()



inputfile = 'preproc_bold-snr-mean' # preproc_bold-edgeview

dat = nib.load(f'/dhcp/scanner_positioning/resliced/{inputfile}.nii.gz')
img = dat.get_fdata()

print(img.shape)
show_slices([img[:,:,z] for z in range(0,dat.shape[2]//2,8)])

plt.savefig(f'/dhcp/scanner_positioning/{inputfile}.jpg')
