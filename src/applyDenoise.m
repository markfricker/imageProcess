function imOut = applyDenoise(imIn, p)
%APPLYDENOISE  Apply deep-learning denoising to a 5-D image array.
%
%   imOut = applyDenoise(imIn, p)
%
% Applies a pre-trained denoising network (currently DnCNN) independently
% to each selected channel/Z/T slice.  Unselected channels are passed
% through unchanged.  If p.use is false the input is returned as-is.
%
% NOTE: named applyDenoise (not denoiseImage) to avoid a name clash with
% MATLAB's built-in denoiseImage(I, net) in the Image Processing Toolbox.
%
% REQUIREMENTS
%   MATLAB Image Processing Toolbox + Deep Learning Toolbox
%
% INPUTS
%   imIn  – [nY nX nC nZ nT] numeric array
%   p     – parameter struct:
%             p.use      – logical; skip if false
%             p.method   – network name string (currently only 'DnCNN')
%             p.channels – [1×nC] logical; which channels to denoise
%
% OUTPUT
%   imOut – denoised array, same size and class as imIn

imOut = imIn;

if ~p.use
    return
end

[~, ~, nC, nZ, nT] = size(imIn);

switch p.method
    case 'DnCNN'
        net = denoisingNetwork('DnCNN');
    otherwise
        error('applyDenoise:unknownMethod', ...
              'Unknown denoising method "%s". Only ''DnCNN'' is supported.', p.method)
end

channels = logical(p.channels(1:nC));

for iT = 1:nT
    for iZ = 1:nZ
        for iC = 1:nC
            if channels(iC)
                im   = double(imIn(:,:,iC,iZ,iT));
                temp = denoiseImage(im, net);   % MATLAB built-in
                imOut(:,:,iC,iZ,iT) = cast(temp, class(imIn));
            end
        end
    end
end
end
