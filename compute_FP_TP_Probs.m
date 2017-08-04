function [FP_Probs, TP_Probs, num_of_tumors, FP_summary, detection_summary] = compute_FP_TP_Probs(results, is_tumor, evaluation_mask, Isolated_Tumor_Cells, EVALUATION_MASK_LEVEL)
% Computes the data for plotting FROC curve and gives a summary of detected
% lesions and false positives
%
% @author: Babak Ehteshami Bejnordi
%
% Input arguments:
% ----------------
% results:                Includes the data collected from the csv file;
% is_tumor:   A boolean variable which is one when the case cotains tumor;
% evaluation_mask:        The labeled evaluation mask;
% Isolated_Tumor_Cells:   An array including the labels containing Isolated Tumor Cells;
% evaluation_mask_level:  The level at which the evaluation mask was made;
%
% Output argument:
% ----------------
% FP_Probs: List of probabilities of false positive findings in the image
%
% TP_Probs: List of probabilities of true positive findings in the image
%
% num_of_tumors: Number of Tumors in the image (excluding Isolate Tumor Cells)
%
% detection_summary:   A matlab structure array with 'fields' that are the labels 
%        of the lesions that should be detected (non-ITC tumors) and 'values'
%        that contain detection details [confidence score, X-coordinate, Y-coordinate]. 
%        Lesions that are missed by the algorithm have an empty value.
%
% FP_summary:   A matlab structure array with 'fields' that represent the 
%        false positive finding number and 'values' that contain detection 
%        details [confidence score, X-coordinate, Y-coordinate]. 
% ----------------
    
    max_label = max(evaluation_mask(:));
    TP_Probs = zeros(1, max_label);
    num_of_tumors = max_label - length(Isolated_Tumor_Cells);
    
    FP_Probs = [];
    detection_summary = struct;
    FP_summary = struct;
   
    counter1 = 1;
    for jj =1:max_label
        if ~any(jj==Isolated_Tumor_Cells)
            label = strcat('Label_ ', int2str(jj));
            detection_summary.(label)= [];
            counter1 = counter1 + 1;
        end
    end

    Counter = 1;    
    if (is_tumor)
        for i = (1: size(results{1},1))
            if (results{3}(i) * results{2}(i) ~= 0)
                HittedLabel = evaluation_mask(results{3}(i)/2^EVALUATION_MASK_LEVEL, results{2}(i)/2^EVALUATION_MASK_LEVEL);
            else
            	HittedLabel = 0;
            end             
            if (HittedLabel == 0)
                FP_Probs(Counter) = results{1}(i); 
                label_num = strcat('FP_ ', int2str(Counter)); 
                FP_summary.(label_num) = [results{1}(i), double(results{2}(i)), double(results{3}(i))];
                Counter = Counter + 1;
            elseif (sum(find(HittedLabel==Isolated_Tumor_Cells))==0)
                %for multiple detections inside one GT object only the one with biggest probability is counted
                if (results{1}(i)>TP_Probs(HittedLabel))
                    TP_Probs(HittedLabel) = results{1}(i);
                    label = strcat('Label_ ', int2str(HittedLabel));
                    detection_summary.(label) = [results{1}(i), double(results{2}(i)), double(results{3}(i))];
                end
            end    
        end 
    else
        for i = (1: size(results{1},1))
            FP_Probs(Counter) = results{1}(i);
            label_num = strcat('FP_ ', int2str(Counter)); 
            FP_summary.(label_num) = [results{1}(i), double(results{2}(i)), double(results{3}(i))];
            Counter = Counter + 1;
        end
    end
end

