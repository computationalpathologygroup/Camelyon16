function [total_FPs, total_sensitivity, FP_summary, detection_summary] = generateFROC(result_dir, masks_dir, EVALUATION_MASK_LEVEL, L0_RESOLUTION)
% The evaluation code for the Camelyon16 challenge on lymph node metastases
% detection
%
% @author: Babak Ehteshami Bejnordi
%
%
% Input:
% -------------------------------------
% result_dir:	Directory containing the CSV files;
% masks_dir:	Directory containing the groundtruth Masks;
% EVALUATION_MASK_LEVEL:  The level at which the evaluation mask is made;
% L0_RESOLUTION:	Pixel resolution of the image at level 0;
%
% Output:
% -------------------------------------
% total_FPs:    an array containing the average number of false positives 
%         per image for different thresholds;
%
% total_sensitivity: an array containing overall sensitivity of the system 
%         for different thresholds;
%
% detection_summary:   A matlab cell array with detection details for each 
%        image. Each row in the cell represent the detection details of a 
%        single image and contains a matlab structure array with 'fields' 
%        that are the labels of the lesions that should be detected 
%        (non-ITC tumors) and 'values' that contain detection details 
%        [confidence score, X-coordinate, Y-coordinate]. Lesions that are 
%        missed by the algorithm have an empty value.
%
% FP_summary:   A matlab cell array which lists the false positive detections  
%        for each image. Each row in the cell represent the detection details 
%        of a single image and contains a matlab structure array with 'fields'   
%        that represent the false positive finding number and 'values' that   
%        contain detection details [confidence score, X-coordinate, Y-coordinate]. 
% -------------------------------------

CSV_file_list = dir(fullfile(result_dir,'*.csv'));
CSV_file_list = {CSV_file_list.name}';

FROC_data = cell(size(CSV_file_list,1),4);
FP_summary = cell(size(CSV_file_list,1),2);
detection_summary = cell(size(CSV_file_list,1),2);

for i=1:numel(CSV_file_list)
    caseID = CSV_file_list{i}(1:end-4);
    fprintf('Evaluating Performance on image: %s\n',caseID);
    is_tumor = strcmp(CSV_file_list{i}(1:5),'Tumor');
    fid = fopen(fullfile(result_dir, CSV_file_list{i}));
    results = textscan(fid,'%f %d32 %d32','Delimiter',',');
    fclose(fid);
        
    if (is_tumor)
        mask_fname = fullfile(masks_dir, strcat(CSV_file_list{i}(1:end-4),'_Mask.tif'));     %# full path to file
        slide = OpenSlide(mask_fname);
        original_mask = slide(:, :, EVALUATION_MASK_LEVEL + 1);
        evaluation_mask = computeEvaluationMask(original_mask, EVALUATION_MASK_LEVEL, L0_RESOLUTION);
        Isolated_Tumor_Cells = computeITCList(evaluation_mask, EVALUATION_MASK_LEVEL, L0_RESOLUTION);
    else
        Isolated_Tumor_Cells = [];
        evaluation_mask = 0;
    end
    
    FROC_data{i, 1} = caseID;
    FP_summary{i, 1} = caseID;
    detection_summary{i, 1} = caseID;
    [FROC_data{i, 2}, FROC_data{i, 3}, FROC_data{i, 4}, FP_summary{i,2}, detection_summary{i,2}] ...
        = compute_FP_TP_Probs(results, is_tumor, evaluation_mask, Isolated_Tumor_Cells, EVALUATION_MASK_LEVEL);  

end

[total_FPs, total_sensitivity] = computeFROC(FROC_data);