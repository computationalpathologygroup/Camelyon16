
% tepisDIR is the directory which holds the MATLAB client for the tEPIS  
% This matlab client should be downloaded from:
% https://github.com/mitkovetta/tepis/tree/master/tepismat
tepisDIR = '...\tepis\tepismat';
addpath(genpath(tepisDIR));

openslide_include_path = '...\openslide\include\openslide';
OpenSlide.initialize(openslide_include_path);

% Directory where all the output .csv files are stored. The csv file names
% should be same as the input image filename.
result_dir = '...\Camelyon16\Results';

% Directory in which the ground truth masks are stored
masks_dir = '...\Ground_Truth\Masks';

EVALUATION_MASK_LEVEL = 5;
L0_RESOLUTION = 0.243;

[total_FPs, total_sensitivity, FP_summary, detection_summary] = generateFROC(result_dir, masks_dir, EVALUATION_MASK_LEVEL, L0_RESOLUTION);
plotFROC(total_FPs, total_sensitivity);