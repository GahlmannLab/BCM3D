% process5Dimages_multiclass process the outputs of niftynet
% thresholding takes custom confidence value and return binary image BW
% numClass denotes which class to process, usually 1 is the background, 2 is  
% interior or surf/sphere ; 3 is boundary or surf+Interior depends on the
% situation; 4 is surf+Interior/rod
% only compatible with nifti_v2
% updated on 20201013 -YW

function [BW]=process5Dimages_multiclass(fileName, thresh, numClass)
    dataFile = fileName;
    image = niftiread_v2(dataFile);
    seg_img = image(:,:,:,1,numClass);
    BW = imbinarize(seg_img,thresh);
end
