function img_intensity_normalization(normalization_parameter)
% fnametif is the tif string, normalization_parameter is the number that
% we multiply to the pixel intensity: to make sure pixel intensities are
% close between each training set.
[fileName,dirName] = uigetfile('*','Load the data');
dataFile = [dirName fileName];
cd (dirName);

try load_nii(fileName);
     imageStack = ans.img * normalization_parameter;
     
catch
     info = imfinfo(fileName);
     warning('try using .tif format');
     imageStack = [];
     numberOfImages = length(info);
     for k = 1:numberOfImages
        currentImage = imread(fnametif, k, 'Info', info);
        currentImage = currentImage.*normalization_parameter;
        flippedImage = fliplr(currentImage);
        Rotateim = imrotate(flippedImage,90);
        imageStack(:,:,k) = Rotateim;
     end
end

datatype = 2;
nii_filled_img1 = make_nii(imageStack , datatype);
newname = [ fileName '.nii'];
save_nii(nii_filled_img1, newname);