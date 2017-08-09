function plotFROC(total_FPs, total_sensitivity)
% This function plots the FROC curve
%
% Input:
% total_FPs:    an array containing the average number of false positives 
%         per image for different thresholds;
%
% total_sensitivity: an array containing overall sensitivity of the system 
%         for different thresholds;

figure, plot(total_FPs, total_sensitivity)
title('Free response receiver operating characteristic curve')
xlabel('Average Number of False Positives')
ylabel('Metastasis detection sensitivity')