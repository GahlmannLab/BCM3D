function [psf_conv, psf_decon, background] = LoadPSF_select(varargin)


%% Load raw PSF
% convolve psf
if varargin{1} == '488psf'
    matObj = matfile('PSFs/PSF488.mat');
    psf_conv = matObj.psf_conv;
    psf_decon = matObj.psf_decon;
    background = matObj.background;
elseif varargin{1} == '561psf'
    matObj = matfile('PSFs/PSF561.mat');
    psf_conv = matObj.psf_conv;
    psf_decon = matObj.psf_decon;
    background = matObj.background;
else
    error('can not read PSF. make sure the PSFs exist and have the correct file names');
end
end