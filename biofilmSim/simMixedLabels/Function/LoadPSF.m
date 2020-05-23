function [psf_rot, psf_decon_rot, background_1] = LoadPSF()
%Process experimental PSF properly, specific for the our home-build 
%lattice light sheet microscope

%% Load raw PSF
% convolve psf
load 488psf.mat psf_raw %Square PSF 20180621
numFrames_1 = size(psf_raw,3);
% deconvolve psf
load 488psf_decon.mat psf_raw_decon %Square PSF 20190221
numFrames_2 = size(psf_raw_decon,3);

%% generate psf, here assume sample step size is 200 nm, psf image acquired by z-step 100nm

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