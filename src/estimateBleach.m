function bleachMean = estimateBleach(imIn, roi, zSource)
%ESTIMATEBLEACH  Measure mean intensity over time to characterise photobleaching.
%
%   bleachMean = estimateBleach(imIn, roi)
%   bleachMean = estimateBleach(imIn, roi, zSource)
%
% Computes the mean intensity within an ROI (or the full frame) for each
% channel and time frame, pooling over Z per `zSource` -- see below. A
% single bleach curve per channel is estimated regardless of nZ: standard
% photobleaching is a global photophysical depletion over TIME, not a
% function of depth, so one curve is the right model here (as opposed to a
% depth-dependent intensity falloff, which is a separate correction).
%
% INPUTS
%   imIn       – [nY nX nC nZ nT] numeric array
%   roi        – [x y w h] ROI, or [] to use the full frame
%   zSource    – (optional) how to pool Z when nZ>1:
%                  'mean'     – average over the whole Z-stack (default).
%                               Pools the most pixels/photons into the
%                               estimate, giving the least noisy curve;
%                               the standard choice for a bleach proxy in
%                               the literature. Can be diluted by largely
%                               out-of-focus/background planes near the
%                               stack edges.
%                  'midplane' – use only the central Z-plane
%                               (round(nZ/2)). Likely to be the most
%                               in-focus, representative single plane, at
%                               the cost of this curve being noisier
%                               (fewer pixels) than 'mean'.
%                  'first'    – use only Z=1 (the old, undocumented
%                               default this function silently had before
%                               2026-10-01 -- kept as an explicit option
%                               for comparison, not recommended).
%                Ignored when nZ==1 (all three are identical then).
%
% OUTPUT
%   bleachMean – [nC × nT] double; mean intensity per channel per frame
%                Returns [] if the ROI is provided but invalid.

if nargin < 3 || isempty(zSource)
    zSource = 'mean';
end

[nY, nX, nC, nZ, nT] = size(imIn);

if isempty(roi)
    % Full-frame measurement
    pos = [1, 1, nX, nY];
else
    pos = roi;
end

if isempty(pos) || ~all(pos(1:4))
    bleachMean = [];
    return
end

x = pos(1); y = pos(2); w = pos(3); h = pos(4);

switch lower(zSource)
    case 'mean'
        zIdx = 1:nZ;
    case 'midplane'
        zIdx = max(1, round(nZ/2));
    case 'first'
        zIdx = 1;
    otherwise
        error('estimateBleach:unknownZSource', ...
            'Unknown zSource "%s". Use ''mean'', ''midplane'' or ''first''.', zSource);
end

roi3D = double(imIn(y:y+h-1, x:x+w-1, :, zIdx, :));  % [h w nC numel(zIdx) nT]

bleachMean = squeeze(mean(roi3D, [1 2 4]));          % pools space + (selected) Z -> [nC × nT]
bleachMean = reshape(bleachMean, nC, nT);
end
