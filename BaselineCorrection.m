% This accompanies "Native LESA-TWIMS-MSI: spatial, conformational and mass
% analysis of proteins and protein complexes"
%%
function [correctedscans] = BaselineCorrection(sumscans,maxb,medb)
% Baseline correction function for synapt data
%   Baseline correction using a moving maximum with a moving median

for i=1:size(sumscans, 2)
    counts = sumscans(:,i);
    
    % Calc the max baseline
    maxbl = maxbaseline(counts,maxb);%25
    
    % Median filter the max baseline to remove the peaks
    % Use a wide filter for this
    medbl = medianbaseline(maxbl, medb);%2000
    
    
    % Subtract the median filter maximum baseline
    counts = counts - medbl;
    counts(counts<0) = 0;
     
    correctedscans(:,i) = counts;
end

end

