function offsets = estimateChannelOffsets(imIn, refChannel, channelMask)
%ESTIMATECHANNEOFFSETS  Estimate translational offsets to align channels.
%
%   offsets = estimateChannelOffsets(imIn, refChannel, channelMask)
%
% Uses intensity-based image registration (monomodal, translation only) on
% a max-Z projection of frame 1 to compute a per-channel XY shift relative
% to the reference channel.
%
% INPUTS
%   imIn        – [nY nX nC nZ nT] numeric array
%   refChannel  – scalar index of the reference (fixed) channel
%   channelMask – [1×nC] logical; true = estimate and apply offset
%
% OUTPUT
%   offsets – struct with fields:
%               offsets.horizontal  – [1×nC] horizontal shift per channel (px)
%               offsets.vertical    – [1×nC] vertical   shift per channel (px)
%               offsets.channels    – channelMask (passed through)
%
% NOTE
%   The returned offsets struct can be passed directly to alignChannels.

nC = size(imIn, 3);
offsets.horizontal = zeros(1, nC);
offsets.vertical   = zeros(1, nC);
offsets.channels   = channelMask;

[optimizer, metric] = imregconfig('monomodal');

% Max-Z projection of frame 1 as the fixed reference
fixed = max(imIn(:,:,refChannel,:,1), [], 4);
tsize = imref2d(size(fixed));
idxC  = find(channelMask);

for k = 1:numel(idxC)
    iC     = idxC(k);
    moving = max(imIn(:,:,iC,:,1), [], 4);
    tform  = imregtform(moving, fixed, 'translation', optimizer, metric);
    offsets.vertical(iC)   = tform.T(3,2);
    offsets.horizontal(iC) = tform.T(3,1);
end
end
