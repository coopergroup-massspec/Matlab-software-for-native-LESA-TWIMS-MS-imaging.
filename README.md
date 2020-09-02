# Matlab-software-for-native-LESA-TWIMS-MS-imaging.

Supplementary information 
Matlab Files
    • DriftScopeSummation.m - script for zero filling and creating the summed counts matrix, as well as extracting the specified region of interest in the drift time dimension.
    • TWIMSIonImage.m - script for processing and creating the ion images, this performs baseline correction and TIC normalisation before creating the image. This requires the SynaptProcessing script to already have been run and requires the following functions
        ◦ BaselineCorrection.m – baseline correction for the spectra
        ◦ maxbaseline.m - moving maximum
        ◦ medianbaseline.m - moving median
        ◦ ProduceFigure.m - Formats ion image
How to run
    1. Download all the files
    2. Change the following inputs in DriftScopeSummation.m
        a. Path – this is the folder location containing the mzML files
        b. imzMLConverterLocation – location of the imzMLConvertor
        c. Scans – number of scans in the mzml file
        d.  Function – what function you want to pull out of the mzml file
        e. SelectionRule – file location of the text file containing the driftscope selection rule
    3. Run SynaptProcessing.m
    4. Change the following variables in TWIMSIonImage.m
        a. maxbwindow – changes size of widow for the moving maximum function
        b. medbwindow – changes size of window for the moving median function
        c. x – number of pixels in the x dimension
        d. y – number of pixels in the y dimension
        e. mzvalue – change this to whatever m/z value you want to show an image of
    5. Run TWIMSIonImage.m
