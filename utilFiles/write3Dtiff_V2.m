function write3Dtiff_V2(array, filename)
% This code modified from B. C. Chen et al. Lattice light-sheet microscopy: 
% Imaging molecules to embryos at high spatiotemporal resolution. 
% Science 346, 1257998 (2014).
outtiff=Tiff(filename, 'w');

dims = size(array);
if length(dims)==2
    dims(3)=1;
end

tagstruct.ImageLength = dims(1);
tagstruct.ImageWidth = dims(2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
if isa(array, 'single')
    tagstruct.BitsPerSample = 32;
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
elseif isa(array, 'uint16')
    tagstruct.BitsPerSample = 16;
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
else
    tagstruct.BitsPerSample = 8;
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
end
tagstruct.SamplesPerPixel = 1;
%tagstruct.RowsPerStrip = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

%outtiff.setTag(tagstruct)

for z=1:dims(3)
    %outtiff.currentDirectory
    outtiff.setTag(tagstruct)
%     outtiff.setTag('ImageWidth', dims(1))
%     outtiff.setTag('ImageLength', dims(2))
%     outtiff.setTag('Photometric',rawtiff.getTag('Photometric'))
%     outtiff.setTag('PlanarConfiguration',rawtiff.getTag('PlanarConfiguration'))
%     outtiff.setTag('BitsPerSample',32)
%     outtiff.setTag('SampleFormat', Tiff.SampleFormat.IEEEFP);
    outtiff.write(array(:,:,z))
    outtiff.writeDirectory()
end
outtiff.close()
