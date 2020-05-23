function [ psf, background] = psf_gen(psf_raw, dz_psf, dz_data, cropToSize, numframes)
% This code modified from B. C. Chen et al. Lattice light-sheet microscopy: 
% Imaging molecules to embryos at high spatiotemporal resolution. 
% Science 346, 1257998 (2014).

%psf_gen resample and crop raw PSF
nz=numframes;
% subtract background estimated from the last Z section
background = mean(mean(psf_raw(:,:,nz)));
psf_raw = double(psf_raw) - background;
% convert all negative pixels to 0
psf_raw(psf_raw<0) = 0.0;

% locate the peak pixel
[~, peakInd] = max(psf_raw(:));
[peakx,peaky,peakz] = ind2sub(size(psf_raw), peakInd);
% crop x-y around the peak
psf_cropped = psf_raw(peakx-cropToSize/2:peakx+cropToSize/2, ...
    peaky-cropToSize/2:peaky+cropToSize/2, :);

% center the PSF in z; otherwise RL decon results are shifted in z
psf_cropped=circshift(psf_cropped, [0, 0, round(nz/2-peakz)]);

[nx,ny,nz]=size(psf_cropped);
% resample PSF to match axial pixel size of raw image
dz_ratio = dz_data / dz_psf;
if dz_ratio > 1
    psf_fft=fftn(psf_cropped);
    new_nz = uint16(round(nz / dz_ratio));
    psf_fft_trunc=complex(zeros(nx,ny,new_nz));
    psf_fft_trunc(:,:,1:new_nz/2)=psf_fft(:,:,1:new_nz/2);
    psf_fft_trunc(:,:,new_nz-new_nz/2+1:new_nz)=psf_fft(:,:,nz-new_nz/2+1:nz); 
    psf=real(ifftn(psf_fft_trunc));
else
    psf = psf_cropped;
end

psf(psf<0) = 0.0;
end

