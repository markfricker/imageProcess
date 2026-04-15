function imOut = filterImageLinear(imIn, p)
%FILTERIMAGELINEAR  Apply a simple linear spatial filter to a 5-D image array.
%
%   imOut = filterImageLinear(imIn, p)
%
% Applies a mean-averaging or Gaussian filter independently to each
% selected channel/Z/T slice.  Designed as a lightweight pre-processing
% step for background estimation; the full adaptive filter pipeline lives
% in DenoiseFilters_sandbox.
%
% INPUTS
%   imIn  – [nY nX nC nZ nT] numeric array
%   p     – parameter struct:
%             p.method       – 'none' | 'mean' | 'Gaussian'
%             p.meanSize     – kernel size for mean filter (odd integer)
%             p.gaussianSize – sigma for Gaussian filter
%             p.channels     – [1×nC] logical; which channels to filter
%
% OUTPUT
%   imOut – filtered array, same size and class as imIn

[nY, nX, nC, nZ, nT] = size(imIn);
imOut    = imIn;
channels = logical(p.channels(1:nC));

if strcmp(p.method, 'none')
    return
end

for iT = 1:nT
    for iZ = 1:nZ
        for iC = 1:nC
            if channels(iC)
                im = imIn(:,:,iC,iZ,iT);
                switch p.method
                    case 'mean'
                        h = fspecial('average', p.meanSize);
                        imOut(1:nY,1:nX,iC,iZ,iT) = imfilter(im, h);
                    case 'Gaussian'
                        imOut(1:nY,1:nX,iC,iZ,iT) = imgaussfilt(im, p.gaussianSize);
                end
            end
        end
    end
end
end
