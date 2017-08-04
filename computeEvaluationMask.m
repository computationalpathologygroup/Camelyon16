function evaluation_mask = computeEvaluationMask(original_mask, EVALUATION_MASK_LEVEL, L0_RESOLUTION)
% This function generates the evaluation masks.
% The evaluation masks are different from the ground truth masks, in that 
% regions in the vicinity of each other are merged together. 
%
% Input arguments:
% ----------------
% original_mask:	Original ground-truth image provided at level of
%                   'evaluation_mask_level';
% evaluation_mask_level:  The level at which the evaluation mask is made;
% L0_resolution:	Pixel resolution of the image at level 0;
%
% Output argument:
% ----------------
% evaluation_mask: The labeled evaluation mask

Threshold = 75/(L0_RESOLUTION * 2 ^ EVALUATION_MASK_LEVEL * 2); %75µm is equivalent to the size of 5 tumor cells
D = bwdist(original_mask(:,:,1));
D = D <= Threshold;
filled_image = imfill(D,'holes');
evaluation_mask = bwlabel(filled_image);
end