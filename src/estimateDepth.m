function depthMean = estimateDepth(imIn, roi)
%ESTIMATEDEPTH  Measure mean intensity per Z-plane to characterise depth-dependent attenuation.
%
%   depthMean = estimateDepth(imIn, roi)
%
% Computes the mean intensity within an ROI (or the full frame) for each
% channel and Z-plane, pooled over the WHOLE time series (all T). Mirrors
% estimateBleach.m with the T and Z roles swapped: depth attenuation
% (absorption/scattering of excitation and emission light with increasing
% depth) is a property of the optical path, not time, so one curve per
% channel is the right model here, the same way one bleach curve per
% channel (not per Z) is right for estimateBleach.
%
% Pooling over all T (rather than a single reference frame) is safe even
% with simultaneous photobleaching: standard bleaching is a uniform
% multiplicative factor across Z (the same assumption correctBleach
% itself uses), so it scales the whole depth curve down evenly without
% distorting its SHAPE -- which is all a depth-attenuation correction
% needs. Depth correction is intended to run BEFORE Bleach correction in
% the pipeline (so Bleach's own Z-pooled estimate is no longer diluted by
% systematic per-depth brightness differences), but estimating it from the
% whole time series rather than a single early frame is correct either way.
%
% INPUTS
%   imIn       – [nY nX nC nZ nT] numeric array
%   roi        – [x y w h] ROI, or [] to use the full frame
%
% OUTPUT
%   depthMean  – [nC × nZ] double; mean intensity per channel per Z-plane
%                Returns [] if the ROI is provided but invalid.
%
% See also: correctDepth, estimateBleach (the time-domain analogue)

[nY, nX, nC, nZ, ~] = size(imIn);

if isempty(roi)
    % Full-frame measurement
    pos = [1, 1, nX, nY];
else
    pos = roi;
end

if isempty(pos) || ~all(pos(1:4))
    depthMean = [];
    return
end

x = pos(1); y = pos(2); w = pos(3); h = pos(4);
roi3D = double(imIn(y:y+h-1, x:x+w-1, :, :, :));  % [h w nC nZ nT]

depthMean = squeeze(mean(roi3D, [1 2 5]));          % pools space + T -> [nC × nZ]
depthMean = reshape(depthMean, nC, nZ);
end
