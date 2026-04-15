function imOut = alignChannels(imIn, offsets)
%ALIGNCHANNELS  Apply pre-computed translational offsets to align channels.
%
%   imOut = alignChannels(imIn, offsets)
%
% Shifts each flagged channel by its stored horizontal and vertical offset
% using an affine2d translation transform.  The output array is the same
% size as the input (pixels shifted out of bounds are lost).
%
% INPUTS
%   imIn    – [nY nX nC nZ nT] numeric array
%   offsets – struct with fields:
%               offsets.horizontal  – [1×nC] horizontal shift per channel (px)
%               offsets.vertical    – [1×nC] vertical   shift per channel (px)
%               offsets.channels    – [1×nC] logical; true = apply shift
%
% OUTPUT
%   imOut   – channel-aligned array, same size and class as imIn

[nY, nX, ~, nZ, nT] = size(imIn);
imOut = imIn;
tsize = imref2d([nY, nX]);
tform = affine2d();
idxC  = find(offsets.channels);

for iT = 1:nT
    for iZ = 1:nZ
        for k = 1:numel(idxC)
            iC = idxC(k);
            tform.T(3,1) = offsets.horizontal(iC);
            tform.T(3,2) = offsets.vertical(iC);
            imOut(:,:,iC,iZ,iT) = imwarp( ...
                imIn(:,:,iC,iZ,iT), tform, ...
                'interp', 'linear', 'OutputView', tsize);
        end
    end
end
end
