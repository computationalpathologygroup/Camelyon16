function Isolated_Tumor_Cells = computeITCList(evaluation_mask, EVALUATION_MASK_LEVEL, L0_RESOLUTION)
% This function defines the labels belonging to isolated tumor cells (ITC)
%
% Description:
%   A region is considered ITC if its longest diameter is below 200µm.
%	As we expanded the annotations by 75µm, the major axis of the object 
%	should be less than 275µm to be considered as ITC (Each pixel is 
%	0.243µm*0.243µm in level 0). Therefore the major axis of the object 
%	in level 5 should be less than 275/(2^5*0.243) = 35.36 pixels.
%
%
% Input arguments:
% ----------------
% evaluation_mask:        The labeled evaluation mask;
% evaluation_mask_level:  The level at which the evaluation mask is made;
% L0_resolution:        Pixel resolution of the image at level 0;
%
% Output argument:
% ----------------
% Isolated_Tumor_Cells: The list of labels belonging to the ITC class
%-----------------

    Threshold = 275/(2 ^ EVALUATION_MASK_LEVEL* L0_RESOLUTION);
    Stats = regionprops(evaluation_mask, 'MajorAxisLength');
    counter = 1;
    Isolated_Tumor_Cells = [];
    for i = 1:size(Stats,1)
        if Stats(i).MajorAxisLength < Threshold
            Isolated_Tumor_Cells(counter) = i;
            counter = counter + 1;
        end
    end 
end