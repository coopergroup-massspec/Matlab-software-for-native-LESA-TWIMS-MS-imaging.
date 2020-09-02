% Copyright 2019 University of Birmingham
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% This file is released as part of the supplementary information for
% O.J. Hale, et al., Native LESA-TWIMS-MSI: spatial, conformational and mass 
% analysis of proteins and protein complexes, XXXX. (20XX),
% DOIXXXXXX
% If you use this software, please acknowledge this work by citing this paper.

% Parts of this code were adapted from SpectralAnalysis, downloaded from 
% https://github.com/AlanRace/SpectralAnalysis on 21/03/2017 and used in this work
% under the conditions of the Apache License, version 2.0.
%% Basline Correction

%change to suit your data
maxbwindow = 25;
medbwindow = 4000;

correctedscans = BaselineCorrection2(ctsmatrix,maxbwindow,medbwindow);

%% Normalisation

% Image size
x = 9;
y = 11;

% Total ion current normalisation
TICfromcorrected = sum(correctedscans,1); % Calculating the total ion count 
NormalisationFactor = (reshape(TICfromcorrected,x,y))'; %reshaping into matrix same size as image

%% Create ion map figure
mzvalue = 2894.7; %MZ values of interes here

for j = 1:size(mzvalue,2) %loops through your mzvalue vector

[f,index] = min(abs(uniquemzs-mzvalue(j))); % finds the mz (and the index) closest to the mzvalue you entered
 
% Summing multiple bins around specified mz (in this case 1 bin eitherm side)
sumpixel = sum(correctedscans(index-1:index+1,:),1);

%Shaping into image dimensions determined by x and y
Ionsum = (reshape(sumpixel,x,y))';

%Creating the normalised image
IonsumNorm = Ionsum./NormalisationFactor;

ProduceFigure(IonsumNorm)
title(strcat('TIC Normalised Baseline Subtract m/z- ', num2str(mzvalue(j))))

end