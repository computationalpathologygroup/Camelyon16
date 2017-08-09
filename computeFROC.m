function [total_FPs, total_sensitivity] = computeFROC(FROC_data)
% This function computes the data needed for plotting FROC curve
%
% Input arguments:
% ----------------
% FROC_data:	A cell array including the List of probabilities of false 
%               positive and true positive findings for all the images; 
%
% Output argument:
% ----------------
% total_FPs:    an array containing the average number of false positives 
%         per image for different thresholds;
%
% total_sensitivity: an array containing overall sensitivity of the system 
%         for different thresholds;
% ----------------

unlisted_FPs = cat(1, [FROC_data{:,2}]);
unlisted_TPs = cat(1, [FROC_data{:,3}]);

all_probs = unique([unlisted_FPs, unlisted_TPs]);
total_FPs = zeros(1, length(all_probs));
total_TPs = zeros(1, length(all_probs));

counter = 1;
for Thresh = all_probs(2:end)
    total_FPs(counter) = sum(unlisted_FPs >= Thresh);
    total_TPs(counter) = sum(unlisted_TPs >= Thresh);
    counter = counter + 1;
end

total_num_of_tumors = sum([FROC_data{:,4}]);

total_FPs = total_FPs./size(FROC_data,1);
total_sensitivity = total_TPs./total_num_of_tumors;