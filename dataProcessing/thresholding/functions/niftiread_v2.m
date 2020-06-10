
%% flip and rotate so that it is consistent with tiff code
function img = niftiread_v2(name)
    temp = niftiread(name);
    img = imrotate(temp,270);
    img = flip(img,2);
end