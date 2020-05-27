# Group-Constrained Subject-Specific (GCSS) Analysis: Make Parcels

Greetings! This is my code to generate a set of parcels from a dataset of individual subject fMRI volumes. Please see the following papers for more on the method:

<li>Fedorenko, E., Hsieh, P.-J., Nieto-Castañon, A., Whitfield-Gabrieli, S. & Kanwisher, N. (2010). A new method for fMRI investigations of language: Defining ROIs functionally in individual subjects. Journal of Neurophysiology, 104(2), 1177-94. DOI: 10.1152/jn.00032.2010.<br> https://pubmed.ncbi.nlm.nih.gov/20410363/</li> 
<li>Julian, J., Fedorenko, E., Webster, J. & Kanwisher, N. (2012). An algorithmic method for functionally defining regions of interest in the ventral visual pathway. Neuroimage, 60(4), 2357-2364. DOI: 10.1016/j.neuroimage.2012.02.055.<br> https://pubmed.ncbi.nlm.nih.gov/22398396/</li>
<br>
And my paper: 
<li>Scott, T.L. and Perrachione, T.K. (2019). Common cortical architectures for phonological working memory identified in individual brains. NeuroImage, 202, 116096. DOI: 10.1016/j.neuroimage.2019.116096.<br> https://pubmed.ncbi.nlm.nih.gov/31415882/</li>

My code is based on the watershed algorithm included in this SPM toolbox developed by Alfonso Nieto-Castañon: https://www.nitrc.org/projects/spm_ss


This code requires Freesurfer and Matlab, and if you want to make pretty pictures, python + a whole bunch of packages (numpy, nibabel, matplotlib, nilearn...). In order to generate a report with the approximate locations of parcels, I use atlases from FSL. Found here -> https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Atlases; license here -> https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Licence


To run: use the Make_My_Parcels.m script to call generate_parcels.m. You will need to specify the full paths to each of your subject maps in a cell array, and a set of options.


I've included part of my data from Scott, T.L. and Perrachione, T.K. (2019) and the parcels generated using that data as an example. These are z-statistic maps of the contrast 4 syllables > 1 syllable during nonword repetition.


I've also included a python script, plot_parcels.py, that will generate a pdf that you can see in the images directory.
