function imOut = registerImage(imIn, tforms)
%REGISTERIMAGE  Apply pre-computed geometric transforms to a 5-D image.
%
%   imOut = registerImage(imIn, tforms)
%
% Applies a per-frame spatial transformation to every channel and Z-plane,
% preserving the original array size.  Transforms are typically estimated
% by a separate registration step (e.g. imregtform on a reference channel).
%
% INPUTS
%   imIn   – [nY nX nC nZ nT] numeric array
%   tforms – nT-element cell array of geometric transformation objects
%            (affine2d, rigid2d, etc.) — one per time frame
%
% OUTPUT
%   imOut  – registered array, same size and class as imIn

[nY, nX, nC, nZ, nT] = size(imIn);
OutputView = imref2d([nY, nX]);
imOut = zeros(size(imIn), 'like', imIn);

for iT = 1:nT
    tform = tforms{iT};
    for iZ = 1:nZ
        for iC = 1:nC
            imOut(:,:,iC,iZ,iT) = imwarp( ...
                squeeze(imIn(:,:,iC,iZ,iT)), tform, ...
                'OutputView', OutputView, 'SmoothEdges', true);
        end
    end
end
end
