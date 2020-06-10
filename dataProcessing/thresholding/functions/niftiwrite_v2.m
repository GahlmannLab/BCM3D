%% flip and rotate so that it is consistent with tiff code
function niftiwrite_v2(img,name)
    newimg = imrotate(img,270);
    newimg = flip(newimg,2);
    niftiwrite(newimg,name);
end