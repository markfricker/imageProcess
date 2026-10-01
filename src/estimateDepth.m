function depthMean = estimateDepth(imIn, roi, mode)
%ESTIMATEDEPTH  Measure per-Z-plane intensity to characterise depth-dependent attenuation.
%
%   depthMean = estimateDepth(imIn, roi)
%   depthMean = estimateDepth(imIn, roi, mode)
%
% Computes an intensity curve across Z, within an ROI (or the full frame),
% for each channel, pooled over the WHOLE time series (all T). Mirrors
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
%   roi        – [x y w h] ROI, or [] to use the full frame. The user is
%                expected to place this ROI over the object of interest --
%                no separate cell/object mask is required or used (none
%                exists yet at this point in the pipeline anyway).
%   mode       – (optional) 'mean' (default) | 'foreground':
%                  'mean'       – average over EVERY pixel in the ROI at
%                                 each Z-plane (original behaviour). Only
%                                 valid when the object's occupancy of the
%                                 ROI is genuinely uniform across Z. For a
%                                 geometrically non-uniform object (e.g. a
%                                 sphere -- more cross-sectional area at
%                                 the equator than the poles), whole-ROI
%                                 averaging conflates true optical
%                                 attenuation with how much of the ROI
%                                 happens to be object at that depth: a
%                                 polar slice's mean gets diluted by
%                                 background pixels for a purely
%                                 geometric reason, not an optical one,
%                                 and "correcting" for that would destroy
%                                 real structure rather than remove a bias.
%                  'foreground' – averages only pixels above a SINGLE
%                                 threshold computed once from the whole
%                                 pooled ROI (all Z, all T) and applied
%                                 IDENTICALLY to every Z-plane. This is
%                                 deliberately NOT a per-Z-plane adaptive
%                                 threshold: an adaptive threshold would
%                                 always select "the brightest fraction of
%                                 whatever is in this plane", which would
%                                 flatten out -- not measure -- the very
%                                 attenuation curve this function exists
%                                 to estimate. Under one fixed threshold, a
%                                 plane with genuinely less/dimmer object
%                                 correctly contributes fewer pixels (or a
%                                 lower mean), removing the occupancy
%                                 confound 'mean' has. Falls back to the
%                                 whole-ROI mean for any Z-plane with zero
%                                 pixels above threshold, rather than NaN.
%                Not a complete fix for every geometry/modality (e.g. a
%                non-optically-sectioned modality where the object's local
%                path-length through the focal slice itself varies with
%                position can still leak through), but removes the
%                dominant occupancy effect.
%
% OUTPUT
%   depthMean  – [nC × nZ] double; intensity per channel per Z-plane
%                Returns [] if the ROI is provided but invalid.
%
% See also: correctDepth, estimateBleach (the time-domain analogue)

if nargin < 3 || isempty(mode)
    mode = 'mean';
end

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
roiVol = double(imIn(y:y+h-1, x:x+w-1, :, :, :));  % [h w nC nZ nT]

switch lower(mode)
    case 'mean'
        depthMean = squeeze(mean(roiVol, [1 2 5]));          % pools space + T -> [nC × nZ]
        depthMean = reshape(depthMean, nC, nZ);

    case 'foreground'
        depthMean = zeros(nC, nZ);
        for iC = 1:nC
            chanVol  = squeeze(roiVol(:,:,iC,:,:));   % [h w nZ nT]
            % ONE threshold from the whole pooled channel volume -- not
            % per-Z -- see header for why a per-plane adaptive threshold
            % would be self-defeating here.
            normVol  = mat2gray(chanVol);
            thresh   = graythresh(normVol);
            for iZ = 1:nZ
                planeNorm = normVol(:,:,iZ,:);
                planeOrig = chanVol(:,:,iZ,:);
                fgMask    = planeNorm > thresh;
                if any(fgMask(:))
                    depthMean(iC,iZ) = mean(planeOrig(fgMask));
                else
                    depthMean(iC,iZ) = mean(planeOrig(:));   % nothing above threshold at this depth
                end
            end
        end

    otherwise
        error('estimateDepth:unknownMode', ...
            'Unknown mode "%s". Use ''mean'' or ''foreground''.', mode);
end
end
