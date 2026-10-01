function [depthMean, depthWeight] = estimateDepth(imIn, roi, mode)
%ESTIMATEDEPTH  Measure per-Z-plane intensity to characterise depth-dependent attenuation.
%
%   depthMean = estimateDepth(imIn, roi)
%   depthMean = estimateDepth(imIn, roi, mode)
%   [depthMean, depthWeight] = estimateDepth(imIn, roi, mode)
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
%   mode       – (optional) 'mean' (default) | 'foreground' |
%                'foreground per section':
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
%                                 IDENTICALLY to every Z-plane. A plane
%                                 with less object contributes fewer
%                                 pixels, removing the occupancy confound
%                                 'mean' has. But as planes dim, only the
%                                 brightest object pixels clear the fixed
%                                 threshold, so for an object with a wide
%                                 intensity range (e.g. ER) the curve
%                                 under-states the dimming (synthetic ER
%                                 stack dimmed to 21%: fit reaches 44%).
%                                 Falls back to the whole-ROI mean for any
%                                 Z-plane with zero pixels above
%                                 threshold, rather than NaN.
%                  'foreground per section'
%                               – Otsu threshold computed separately for
%                                 each Z-plane (pooled over T), averaging
%                                 the RAW intensities above it. The
%                                 threshold only selects pixels; their raw
%                                 mean still scales with the plane's
%                                 dimming, so the curve is not flattened,
%                                 and the same share of the object is kept
%                                 at every depth (same ER stack: 24% vs a
%                                 true 21%). A plane past the end of the
%                                 object is pure noise, which Otsu still
%                                 splits in two: a plane whose foreground-
%                                 minus-background mean is under 5% of the
%                                 largest in the stack is treated as empty
%                                 (weight eps, like 'foreground's
%                                 fallback).
%                Not a complete fix for every geometry/modality (e.g. a
%                non-optically-sectioned modality where the object's local
%                path-length through the focal slice itself varies with
%                position can still leak through), but removes the
%                dominant occupancy effect.
%
% OUTPUT
%   depthMean   – [nC × nZ] double; intensity per channel per Z-plane
%                 Returns [] if the ROI is provided but invalid.
%   depthWeight – [nC × nZ] double; number of pixels behind each
%                 depthMean entry (w*h*nT for 'mean', always the same at
%                 every Z; the foreground pixel count for the foreground
%                 modes,
%                 which can vary hugely by Z-plane -- pass this to
%                 correctDepth so its exponential fits weight a precise,
%                 many-pixel estimate more than a noisy few-pixel one,
%                 rather than treating every Z-plane's mean as equally
%                 trustworthy regardless of how many pixels it came from).
%
% See also: correctDepth, estimateBleach (the time-domain analogue)

if nargin < 3 || isempty(mode)
    mode = 'mean';
end

[nY, nX, nC, nZ, nT] = size(imIn);

if isempty(roi)
    % Full-frame measurement
    pos = [1, 1, nX, nY];
else
    pos = roi;
end

if isempty(pos) || ~all(pos(1:4))
    depthMean   = [];
    depthWeight = [];
    return
end

x = pos(1); y = pos(2); w = pos(3); h = pos(4);
roiVol = double(imIn(y:y+h-1, x:x+w-1, :, :, :));  % [h w nC nZ nT]

switch lower(mode)
    case 'mean'
        depthMean = squeeze(mean(roiVol, [1 2 5]));          % pools space + T -> [nC × nZ]
        depthMean = reshape(depthMean, nC, nZ);
        depthWeight = repmat(w*h*nT, nC, nZ);                % every Z-plane uses the same pixel count

    case 'foreground'
        depthMean   = zeros(nC, nZ);
        depthWeight = zeros(nC, nZ);
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
                nFg       = nnz(fgMask);
                if nFg > 0
                    depthMean(iC,iZ)   = mean(planeOrig(fgMask));
                    depthWeight(iC,iZ) = nFg;
                else
                    % Nothing above threshold at this depth: the fallback
                    % mean is pure background/noise, not a measurement of
                    % the attenuation curve's value here -- give it a
                    % near-zero weight (not exactly zero, to avoid an
                    % all-zero weight vector if this happens at every Z)
                    % so a weighted fit effectively ignores it, rather
                    % than the opposite mistake of treating "the whole
                    % plane went into this mean" as if it were the MOST
                    % trustworthy point.
                    depthMean(iC,iZ)   = mean(planeOrig(:));
                    depthWeight(iC,iZ) = eps;
                end
            end
        end

    case 'foreground per section'
        depthMean   = zeros(nC, nZ);
        depthWeight = zeros(nC, nZ);
        minContrast = 0.05;   % of the stack's largest fg-bg difference; below = no object
        for iC = 1:nC
            contrast = zeros(1, nZ);
            for iZ = 1:nZ
                planeOrig = double(squeeze(roiVol(:,:,iC,iZ,:)));   % [h w nT]
                fgMask    = planeOrig > graythresh(mat2gray(planeOrig)) * ...
                    (max(planeOrig(:)) - min(planeOrig(:))) + min(planeOrig(:));
                if any(fgMask(:)) && ~all(fgMask(:))
                    depthMean(iC,iZ)   = mean(planeOrig(fgMask));
                    depthWeight(iC,iZ) = nnz(fgMask);
                    contrast(iZ)       = depthMean(iC,iZ) - mean(planeOrig(~fgMask));
                else
                    depthMean(iC,iZ)   = mean(planeOrig(:));   % flat plane: nothing to split
                end
            end
            % A plane of pure noise is still split in two by Otsu, but its
            % fg-bg difference is a small fraction of an object plane's.
            isEmpty = contrast < minContrast * max(contrast);
            depthWeight(iC, isEmpty) = eps;
        end

    otherwise
        error('estimateDepth:unknownMode', ...
            'Unknown mode "%s". Use ''mean'', ''foreground'' or ''foreground per section''.', mode);
end
end
