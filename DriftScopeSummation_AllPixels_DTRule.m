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

%% Inputs

Path='C:\Users\haleoj\Documents\Data\Synapt G2-S\June 2019\2019_06_24_TWIMSimage\'; %folder containing the mzML files

imzMLConverterLocation='C:\Users\haleoj\OneDrive\NAMS\imzMLConverter_1.3\imzMLConverter.jar\' ; %location of the imzml converter

filestoprocess=dir([Path '*.mzML']); %creates a list of the mzML files in folder 'Path'

%FIND 
Scans = 259; 

Function = 2; % which function you want to pull out from the mzml file

%allows the use of the imzml converter
javaclasspath(imzMLConverterLocation);
javaaddpath(imzMLConverterLocation);


%% Picking out the right function (1 or 2)
%only performing on one pixel, Rian's pixels are all the same size file so
%the indexing is the same, if not you will have to run this section for
%evey pixel.

index = []; %creating an empty matrix in which to put the indexes of the spectrums that are of a specific function
counter = 1;

%change, this is the pixel number
filename = filestoprocess(1).name;  %chose p26 to do this on
filePath = [Path filesep filename];
mzML = imzMLConverter.MzMLHandler.parsemzML(filePath);

for j = 0:Scans-1 % for loop running through all the individual scans in the file
    
    spectrum = mzML.getRun().getSpectrumList().getSpectrum(j); %pulls out spectum j
    ID = char(spectrum.getID()); % pulls out the spectrum ID from the mzml
    ID = str2double(ID(10)); % the 10th character of spectrum ID (function 1 or 2)
    
    %if statement - if function equals chosen value pull out the index of
    %the spectrum
    if ID == Function %choosing the function (input values 1 or 2)
        index(counter) = j; %putting ID into index matrix
        counter = counter + 1;
        
    else
        
    end
    
end

%% MZ list from the scans with drift time data

totalmzs = []; %empty matix to put the total mzs in

for z=1:size(filestoprocess,1)  %looping through the pixel files
    
    
    filename = filestoprocess(z).name; %reads the name of the working file from the list 'filestoprocess'
    
    filePath = [Path filesep filename]; %generates the individual file path by adding the working file name to the working directory path
    
    mzML = imzMLConverter.MzMLHandler.parsemzML(filePath);
    
    mzArray = cell(size(index,2),1); %cell(Scansend(i)-Scansstart(i)+1,1);
    
    counter2 = 1;
    
    
    for k = index   %loop through the scans in a pixel
        
        spectrum = mzML.getRun().getSpectrumList().getSpectrum(k); %calling out the spectrum from the scan specified
        
        if isempty(spectrum)
            
            disp('Warning: Spectrum is empty :(')
        else
            mzArray{counter2} = spectrum.getmzArray(); %gettin g the mz array
        end
        counter2 = counter2+1; %moving the counter along one
        
    end %scans in pixels
    
    total = cat(1,mzArray{:}); %concatonating the mz arrays for each scan of the pixel
    totalmzs = cat(1,totalmzs,total); %concatonating the mz arrays from the pixels with the previous mz arrays
    
end %pixel files

uniquemzs = unique(totalmzs); %finding the unique mz values in the total mz, this will be needed later


%% Counts
%change
SelectionRule = dlmread('C:\Users\haleoj\OneDrive\NAMS\Driftscope Selection Rules\ROI_proteins_kidney_f43simp.txt','',6,0); %Pulls in the rule (text file)

%Counts matrix
ctsmatrix_DT=zeros(size(uniquemzs,1),size(filestoprocess,1)); % creating an empty counts matrix to store the counts in (mz size x 99 pixels)
i=1; %initialising i will be used to cycle through the 200 drift bins

for a=1:size(filestoprocess,1)
    
    tempcts = zeros(size(uniquemzs,1),200);
    
    filename = filestoprocess(a).name; %reads the name of the working file from the list 'filestoprocess'
    
    filePath = [Path filesep filename]; %generates the individual file path by adding the working file name to the working directory path
    
    mzML = imzMLConverter.MzMLHandler.parsemzML(filePath);
    
    
    for z = index %pulling out the spectrums from function that was specified earlier in one pixel
        
        spectrum = mzML.getRun().getSpectrumList().getSpectrum(z); %call the spectrum out
        if isempty(spectrum)
        disp('Warning: Spectrum is still empty :(')
        else
        mzArray = spectrum.getmzArray(); %pull out the mz array
        countz = spectrum.getIntensityArray(); %pull the respective counts array
        mzindex = ismember(uniquemzs,mzArray); % pi;;
        end
        s = 1;
        
        for n=1:length(uniquemzs) %looping through the length of the unique mz list
            
            if mzindex(n) == 1 % if the nth element of the logic matrix=1
                tempcts(n,i) = countz(s)+tempcts(n,i); %put the associating count in the nth element of the counts matrix
                s = s+1;
            else %else do nothing, counter stays the same
            end %end of if statement
            
        end %end of loop through the mz list
        
        
        if i == 200
            i  = 1;
        else
            i = i+1;
        end
        
    end %end of loop through index
    
    
    %Extraction of region of interest
    driftcounts = zeros(size(tempcts)); %creates a matrix of zeros the size of the counts matrix
    
    for d = 1:200 % looping through the drift bins
        
        MzFrom = SelectionRule(d,2); % Pulling out the second colum of the text file which corresponds to the mz start value
        
        if MzFrom ~= -1 % if the mz from value is not -1 do the following
            
            MzTo = SelectionRule(d,3); % pull out end mz value for the correspondinf drift bin
            
            [f, CountsStart] = min(abs(uniquemzs-MzFrom)); % finding the closest value to the mz start value from the overall mz matrix
            [g, CountsEnd] = min(abs(uniquemzs-MzTo)); % finding the closest value to the mz end value from the overall mz matrix
            
            driftcounts(CountsStart:CountsEnd,d) =  tempcts(CountsStart:CountsEnd,d); %pulling out the part of the patrix that xoressponds to this drift time and mzvalues and putting into driftcounts matrix
            
        else
            
        end %end of if statement
        
    end %end of for loop
    
    ctsmatrix_DT(:,a) = sum(driftcounts,2);
    
end %end of loop through pixel files


% ctsmatrix will be required later







