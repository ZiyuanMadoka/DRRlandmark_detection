%% read and edit DICOM RP files
cd('H:/Desktop/diagnositc/AMC_53/20170921Diag2013Feb/')
fileset = dir('H:/Desktop/diagnositc/AMC_53/20170921Diag2013Feb/RP*');
X=dicomread(fileset(2).name);