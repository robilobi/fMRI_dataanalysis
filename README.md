ANALYSIS of fMRI data with statistical parametric mapping (SPM12; Welcome Trust Centre for Neuroimaging; http:// www.fil.ion.ucl.ac.uk/spm/software/spm12) using standard spatial pre-processing procedures. 

PREPROCESSING (preprocessing_firstlevel.m)
- slice time correction (using cubic spline interpolation),
- spatial realignment, 
- co-registration of functional and anatomical data, 
- spatial normalization into the MNI (Montreal Neurological Institute) stereotactic space, that includes resampling to 2 × 2 × 2 mm voxel size. 
-  spatial low-pass filtering using a 3D Gaussian kernel with full-width at half-maximum (FWHM) of 8 mm, 
- temporal high-pass filtering with a cut-off of 1/128 Hz to eliminate low-frequency drifts. 

FIRSTLEVEL STATISTICS (preprocessing_firstlevel.m)
- set up t-contrasts of interest at the first level
- the evoked hemodynamic response to the onset of event of interest is modeled for each condition as boxcars convolved with a hemodynamic response function (HRF). 
- estimated motion re-alignment parameters are added to the design as covariates of no interest to regress out residual motion artifacts and increase statistical sensitivity.

SECONDLEVEL STATISTICS (secondlevel_ANOVAfullfactorial.m)


!! For statistical thresholding, a Monte Carlo simulation implemented in Matlab (1000 iterations, no volume mask) can be run to define the minimum cluster extent threshold of  a voxel-level uncorrected p-value of .001 to ensure whole-volume type I error probability smaller than .05 (Slotnick et al., 2003; code available at http://www2.bc.edu/sd-slotnick/scripts.html). 
