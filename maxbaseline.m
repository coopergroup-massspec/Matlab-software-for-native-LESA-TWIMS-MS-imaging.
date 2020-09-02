% This accompanies "Native LESA-TWIMS-MSI: spatial, conformational and mass
% analysis of proteins and protein complexes"

function [ baseline ] = maxbaseline( counts, hw )
%MAXBASELINE 

    baseline = zeros(size(counts));
    
    for i=1+hw:length(counts)-hw
        baseline(i) = max(counts(i-hw:i+hw));
    end
    
    for i=1:hw
        baseline(i) = max(counts(1:i+hw));
    end
    
    for i=length(counts)-hw+1:length(counts)
        baseline(i) = max(counts(i-hw:length(counts)));
    end
    
end

