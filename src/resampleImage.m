function imOut = resampleImage(imIn, p)
%RESAMPLEIMAGE  Spatially downsample (or upsample) a 5-D image array.
%
%   imOut = resampleImage(imIn, p)
%
% Downsamples by p.factor using bilinear interpolation.  For downsampling
% (factor > 1) a guided-filter anti-aliasing pass is applied first to each
% slice.  Binary (logical) images use nearest-neighbour interpolation and
% skip the anti-alias step.
%
% INPUTS
%   imIn  – [nY nX nC nZ nT] numeric or logical array
%   p     – parameter struct:
%             p.factor – scalar resampling factor (>1 = downsample)
%
% OUTPUT
%   imOut – resampled array (same class as imIn)

factor = double(p.factor);

if islogical(imIn)
    imOut = imresize(imIn, 1/factor, 'nearest');
    return
end

[~, ~, nC, nZ, nT] = size(imIn);

% Anti-alias kernel: odd, minimum 3
kernelSize = max(3, 2*floor(3*factor/2) + 1);

% Pre-compute output dimensions from a single test slice
testSlice = imresize(imIn(:,:,1,1,1), 1/factor, 'bilinear');
[nYr, nXr] = size(testSlice, [1 2]);

imOut = zeros(nYr, nXr, nC, nZ, nT, 'like', imIn);

for iT = 1:nT
    for iZ = 1:nZ
        for iC = 1:nC
            slice = imIn(:,:,iC,iZ,iT);
            if factor > 1
                slice = imguidedfilter(slice, slice, ...
                    'NeighborhoodSize', [kernelSize kernelSize]);
            end
            imOut(:,:,iC,iZ,iT) = imresize(slice, 1/factor, 'bilinear');
        end
    end
end
end
