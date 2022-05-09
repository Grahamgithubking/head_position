# Converts six FSL motion parameters to FLIRT format transformation matrix
#  Provide filename of motion corrected series, used to adjust origin for correct rotation
#  Rhodri Cusack, TCIN Dublin, 2021-11-24, cusackrh@tcd.ie
import numpy as np
import sys
import nibabel as nib
import argparse
import math

parser = argparse.ArgumentParser(description='Convert motion parameters to FLIRT transformation matrix')
parser.add_argument('nifti',  type=str, nargs=1, help='nifti file to define space of transformation')
parser.add_argument('movepar',  type=float, nargs=6, help='six motion parameters (translations first)')
parser.add_argument('--invert', action='store_true', help='invert final transformation matrix')
parser.add_argument('--degrees',  action='store_true', help='rotations in degrees not radians')

args = parser.parse_args()

# Centre of rotation for MCFLIRT is centre of image, not voxel (0,0,0) as in FSL definition of matrix
#  https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT/FAQ#What_is_the_format_of_the_matrix_used_by_FLIRT.2C_and_how_does_it_relate_to_the_transformation_parameters.3F
# Load up volume
img = nib.load(args.nifti[0])
af = img.affine
origin_mm = img.affine[:4, 3]
com_vox = np.ones((4))
com_vox[:3] = img.header['dim'][1:4]/2  # centre of mass in voxels
com_mm = img.affine@com_vox.T

# Convert arguments to double
p = np.array([np.double(x) for x in args.movepar])

if args.degrees:
    p[3:] = math.pi* p[3:] / 180


# Set up transformation matrices
originoff = np.eye(4)   # adjust origin to centre of mass
originoff[:3, 3] = com_mm[:3] - origin_mm[:3]

# To mimic MCFLIRT values, need to flip  x
originoff[0,:] = -originoff[0,:]

# Create translations and rotations
t = np.array([[1, 0, 0, p[0]], [0, 1, 0, p[1]], [0, 0, 1, p[2]], [0, 0, 0, 1]])
r1 = np.array([[1, 0, 0, 0], [0, np.cos(p[3]), np.sin(p[3]), 0],
               [0, -np.sin(p[3]), np.cos(p[3]), 0], [0, 0, 0, 1]])
r2 = [[np.cos(p[4]), 0, np.sin(p[4]), 0], [0, 1, 0, 0],
      [-np.sin(p[4]), 0, np.cos(p[4]), 0], [0, 0, 0, 1]]
r3 = [[np.cos(p[5]), -np.sin(p[5]), 0, 0], [np.sin(p[5]),
                                           np.cos(p[5]), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
originoffinv = np.linalg.inv(originoff) # adjust origin back
rot = r1@r2@r3

# Creatre final matrix
#  Translation before flip of x axis
mat = t@originoff@rot@originoffinv

# Return inverse
if args.invert:
    mat=np.linalg.inv(mat)

for row in mat:
    for item in row:
        print('%0.6f ' % item, end='')
    print('')
