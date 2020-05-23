function [psf_rot, psf_decon_rot, background_1] = LoadPSF_select(varargin)


%% Load raw PSF
% convolve psf
if varargin{1} == '488psf'
    % raw psf
    matObj = matfile('PSFs/488psf.mat'); %Square PSF 20180621
    details = whos(matObj);
    psfName = details.name;
    psf_raw = matObj.(psfName);
    numFrames_1 = size(psf_raw,3);
    % deconvolved psf
    matObj = matfile('PSFs/488psf_decon.mat'); %Square PSF 20180621
    details = whos(matObj);
    psfName = details.name;
    psf_raw_decon = matObj.(psfName); %Square PSF 20190221
    numFrames_2 = size(psf_raw_decon,3);
elseif varargin{1} == '561psf'
    % raw psf
    matObj = matfile('PSFs/561psf.mat'); %Square PSF 20180621
    details = whos(matObj);
    psfName = details.name;
    psf_raw = matObj.(psfName); %561 Square PSF 20190618
    numFrames_1 = size(psf_raw,3);
    % deconvolved psf
    matObj = matfile('PSFs/561psf_decon.mat'); %Square PSF 20180621
    details = whos(matObj);
    psfName = details.name;
    psf_raw_decon = matObj.(psfName); %Square PSF 20190221
    numFrames_2 = size(psf_raw_decon,3);
else
    error('can not read PSF. make sure the PSFs exist and have the correct file names');
end
%% generate psf, here assume sample step size is 200 nm, psf image acquired
% by z-step 100nm
rotateByAngle=31.8;
dz_data_ratio = sin(rotateByAngle*pi/180);
dz_step = 200;
dz_psf = 100;
dz_data = dz_step*dz_data_ratio;
cropToSize = 48;
[psf, background_1] = psf_gen(psf_raw, dz_psf, dz_data, cropToSize, numFrames_1);
%background_1 get from psf for convolution
zx_aspratio = 1;
psf_rot = rotate3D( psf, rotateByAngle, zx_aspratio);
[psf_decon, ~] = psf_gen(psf_raw_decon, dz_psf, dz_data, cropToSize, numFrames_2);
psf_decon_rot = rotate3D( psf_decon, rotateByAngle, zx_aspratio);
end