% This accompanies "Native LESA-TWIMS-MSI: spatial, conformational and mass
% analysis of proteins and protein complexes"

function [ baseline ] = medianbaseline( counts, hw )
%MEADIANBASELINE 

    baseline = zeros(size(counts));
    
    for i=1+hw:length(counts)-hw
        baseline(i) = median(counts(i-hw:i+hw));
    end
    
    for i=1:hw
        baseline(i) = median(counts(1:i+hw));
    end
    
    for i=length(counts)-hw+1:length(counts)
        baseline(i) = median(counts(i-hw:length(counts)));
    end
    
end

