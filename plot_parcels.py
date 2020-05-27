#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 10:25:25 2020

@author: Terri L. Scott (tlscott@bu.edu)

Creates a vector image of your parcels on a glass brain. I plot them just as 
outlines, because they are not plotted according to what should actually be in 
front of what, so you need to fix that by hand for a publication-ready figure. 
I import the resulting pdf from this script into inkscape and customize from
there.

"""

from nilearn import plotting
from nilearn import image
import numpy as np
import matplotlib
from matplotlib import cm
import nibabel as nib

my_colormap = 'Spectral'
PROJECT_DIR = '/Users/tlscott/make_parcels/'

# Path to your volume of parcels
PARCEL_DIR = PROJECT_DIR + 'parcels/My_Experiment_26-May-2020/'
# Path to where you want to save the image
IMG_PATH = PROJECT_DIR + 'images/'
IMG_NAME = 'parcels.pdf'
# Full path to parcel volume with filename
PARCEL_VOL = PARCEL_DIR + 'My_Experiment_probability_map_thresh2subjs_smoothed_parcels_sig.nii.gz'

# Loads the volume as a nifti object
img = image.index_img(PARCEL_VOL,0)
# Gets just the volume data, converting the volume to an easy to manipulate 3D array
img_array = img.get_data() # This line will need to be updated for future versions of nilearn
# Each parcel has a unique number. 
img_vals = np.unique(img_array)

# Make a custom color map with as many values as there are parcels
a_colormap = cm.get_cmap(my_colormap,len(img_vals)) 

# Opens the display
display = plotting.plot_glass_brain(None, display_mode='lzr')

for i in img_vals:
    # For each unique value not 0, plot contour
    if i > 0:
        # Get only one parcel at a time
        temp_img = np.copy(img_array)
        temp_img[temp_img != i] = 0
        # Convert the array back to a nifti object
        temp_img_nii = nib.Nifti1Image(temp_img[:,:,:],img.affine)
        # Assign a color (can be any value btwn 0 and 1)
        temp_color = a_colormap(i/img_vals[-1])
        # Color format conversion nonsense. Nilearn only likes hex values.
        my_color = matplotlib.colors.rgb2hex(temp_color[0:3])  
        # Finally plot the contour
        display.add_contours(temp_img_nii, levels=[0.5], colors=my_color, linewidths=3.)
            
# Don't forget to close the display
display.savefig('%s/%s' %(IMG_PATH,IMG_NAME))        
display.close() 

