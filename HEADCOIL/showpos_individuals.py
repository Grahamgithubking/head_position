import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt

def show_slices(slices):
    """ Function to display row of image slices """




nsub=10
nslices=12

inputfile = 'preproc_bold-snr-merged' # 'preproc_bold-snr-mean' preproc_bold-edgeview
dat = nib.load(f'/dhcp/scanner_positioning/resliced/{inputfile}.nii.gz')
img = dat.get_fdata()

fig, axes = plt.subplots(nsub, nslices)

x0=1
x1=128
y0=1
y1=128

for subind in range(nsub):
    slices=[img[x0:x1,y0:y1,z,subind] for z in range(0,dat.shape[2]//2,8)]

    for i, slice in enumerate(slices):
        axes[subind][i].imshow(slice.T, cmap="gray", origin="lower", vmin=0, vmax=100)
        axes[subind][i].set_axis_off()

    plt.savefig(f'/dhcp/scanner_positioning/{inputfile}.jpg')