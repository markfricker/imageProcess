function [imOut, results] = correctDepth(imIn, depthMean, p, depthWeight)
%CORRECTDEPTH  Correct depth-dependent intensity attenuation from a Z-stack image.
%
%   [imOut, results] = correctDepth(imIn, depthMean, p)
%   [imOut, results] = correctDepth(imIn, depthMean, p, depthWeight)
%
% Estimates per-plane, per-channel correction factors from a depth-
% attenuation curve (depthMean from estimateDepth) and applies them to the
% image. Output intensities are normalised to [0,1] via mat2gray before
% scaling. Mirrors correctBleach.m with the T and Z roles swapped.
%
% ONE correction factor per (channel, Z-plane) is applied UNIFORMLY to
% every T and every pixel in that plane: depth attenuation is a static
% property of the optical path (absorption/scattering increasing with
% depth), not a function of time, so a single curve per channel is the
% right model (the mirror image of correctBleach applying one per-(channel,
% frame) factor uniformly across every Z-plane).
%
% METHODS (identical menu to correctBleach)
%   'ratio'            – correction = 1 / normalised_mean(z)
%   'exponential'      – fit single exponential (exp1), correct by 1/fit(z)
%   'double exponential'– fit double exponential (exp2), correct by 1/fit(z)
%   'histogram match'  – match each plane's histogram to Z=1 (no fitting)
%
% REQUIREMENTS
%   'exponential' and 'double exponential' require Curve Fitting Toolbox.
%
% INPUTS
%   imIn        – [nY nX nC nZ nT] numeric array
%   depthMean   – [nC × nZ] double from estimateDepth (one curve per
%                 channel, pooled over the whole time series)
%   p           – parameter struct:
%                   p.method     – method string (see above)
%                   p.norm2Start – logical; true = normalise to Z=1
%                                  (shallowest plane), false = normalise to
%                                  the last (deepest) plane
%                   p.channels   – [1×nC] logical; which channels to correct
%   depthWeight – (optional) [nC × nZ] double, the depthWeight output of
%                 estimateDepth (pixel count behind each depthMean entry).
%                 Used as 'exponential'/'double exponential' fit weights,
%                 since Var(mean) ~ 1/n -- a Z-plane whose mean came from
%                 many pixels (e.g. the equator of a sphere under
%                 'foreground' mode) is a more precise estimate than one
%                 from a handful of pixels (near the poles) and should
%                 pull the fitted curve harder. Omit, or pass [], for
%                 uniform weighting (every Z-plane trusted equally) --
%                 always correct for 'mean' mode, where every plane uses
%                 the same pixel count anyway. Not used by 'ratio' or
%                 'histogram match', neither of which fits across points.
%
% OUTPUTS
%   imOut   – corrected array, same size as imIn, class determined by
%             original class (uint8 / uint16 / single)
%   results – struct with fields:
%               results.correctionFactor – [nC × nZ] applied per-plane factors
%               results.norm             – [nC × nZ] normalised depth curve
%               results.fitObjects       – {1×nC} cfit objects (exp methods only)
%               results.gof              – {1×nC} goodness-of-fit (exp methods only)
%               results.normRefWarning   – '' if fine, else a message naming
%                                          which channel(s) have an unstable
%                                          normalisation reference (see below)
%
% NORMALISATION REFERENCE STABILITY (p.norm2Start)
%   'ratio'/'exponential'/'double exponential' all normalise the curve by
%   dividing by the reference plane's value (Z=1 or Z=end) before fitting/
%   inverting. If that reference plane happens to be close to zero (e.g.
%   genuinely no signal there -- an out-of-focus or background-only edge
%   plane), every OTHER plane's correctionFactor = normRef/value collapses
%   toward zero, crushing the corrected image everywhere except at the
%   reference plane itself -- not a subtle effect, a "daft multiplier".
%   This is detected (reference < 10% of the channel's own peak) and
%   returned in results.normRefWarning (and raised via warning()) rather
%   than silently producing a badly-scaled correction; it does not block
%   the correction (the caller decides what to do -- e.g. offer to
%   normalise to the other end instead).
%
% See also: estimateDepth, correctBleach (the time-domain analogue)

[~, ~, nC, nZ, ~] = size(imIn);
channels = logical(p.channels(1:nC));
z        = (1:nZ)';
origClass = class(imIn);

if nargin < 4 || isempty(depthWeight)
    depthWeight = ones(nC, nZ);
end

% Single-plane (non-Z-stack) input: depth correction is a true no-op here,
% not just "nothing to correct" -- the fit-based methods need at least a
% couple of Z-planes to fit a decay curve to and would error on one point.
if nZ <= 1
    imOut = imIn;
    results.correctionFactor = ones(nC, 1);
    results.norm             = ones(nC, 1);
    results.fitObjects       = cell(1, nC);
    results.gof              = cell(1, nC);
    results.normRefWarning   = '';
    return
end

% Normalise depth curve
if p.norm2Start
    normRef = depthMean(:, 1);
    refLabel = 'first (Z=1)';
else
    normRef = depthMean(:, end);
    refLabel = 'last (deepest)';
end
normRefRaw = normRef;   % pre-guard, for the stability check below
normRef(normRef == 0) = 1;   % guard against exact-zero reference (division)
norm = depthMean ./ normRef;

% See header: a reference point close to zero (not just exactly zero)
% still makes every other plane's correction factor collapse toward zero
% once inverted. Flag it rather than silently returning a badly-scaled
% correction -- it does not block the method the caller asked for.
normRefWarning = '';
if ~strcmp(p.method, 'histogram match')
    peakVal = max(depthMean, [], 2);
    smallRefFrac = 0.1;   % reference point is < 10% of this channel's own peak
    isBadRef = channels(:) & (normRefRaw < smallRefFrac * peakVal);
    if any(isBadRef)
        normRefWarning = sprintf(['Depth correction: channel(s) %s -- the %s-plane reference value is ' ...
            'less than %.0f%% of that channel''s peak intensity. Normalising to it will make the ' ...
            'correction factor collapse toward zero everywhere else, crushing the corrected image. ' ...
            'This usually means that plane has little/no real signal (e.g. an out-of-focus or ' ...
            'background-only edge plane) -- this correction is probably not appropriate for this ' ...
            'curve as configured; consider normalising to the other end instead.'], ...
            mat2str(find(isBadRef)'), refLabel, 100*smallRefFrac);
        warning('correctDepth:unstableNormReference', '%s', normRefWarning);
    end
end

correctionFactor = ones(nC, nZ);
fitObjects       = cell(1, nC);
gof              = cell(1, nC);

% Working copy as double for mat2gray. mat2gray's min/max is taken per
% (channel, T-frame) across the WHOLE Z-stack, so the planes keep their
% relative brightness -- the depth curve the factors below correct. (The
% mirror of correctBleach, which scales per (channel, Z-plane) across T.)
% Scaling each Z-plane separately would already flatten most of the
% attenuation, and the factors would then over-correct the deep planes.
nT = size(imIn, 5);
temp = zeros(size(imIn));
for iT = 1:nT
    for iC = 1:nC
        temp(:,:,iC,:,iT) = mat2gray(imIn(:,:,iC,:,iT));
    end
end

switch p.method
    case 'ratio'
        for iC = 1:nC
            correctionFactor(iC,:) = 1 ./ norm(iC,:);
        end

    case 'exponential'
        for iC = 1:nC
            [fitObjects{iC}, gof{iC}] = fit(z, norm(iC,:)', 'exp1', 'Weights', depthWeight(iC,:)');
            correctionFactor(iC,:)    = 1 ./ feval(fitObjects{iC}, z)';
        end

    case 'double exponential'
        for iC = 1:nC
            [fitObjects{iC}, gof{iC}] = fit(z, norm(iC,:)', 'exp2', 'Weights', depthWeight(iC,:)');
            correctionFactor(iC,:)    = 1 ./ feval(fitObjects{iC}, z)';
        end

    case 'histogram match'
        % Each Z-plane is matched to the Z=1 reference (pooled over T),
        % not a single (Z,T) frame.
        ref = squeeze(mean(double(imIn(:,:,:,1,:)), 5));   % [nY nX nC]
        for iZ = 1:nZ
            for iC = 1:nC
                if channels(iC)
                    for iT = 1:nT
                        temp(:,:,iC,iZ,iT) = mat2gray( ...
                            imhistmatch(imIn(:,:,iC,iZ,iT), ref(:,:,iC)));
                    end
                end
            end
        end
        tempMean = squeeze(mean(temp, [1 2 5]));    % pooled over space+T -> [nC × nZ]
        tempMean = reshape(tempMean, nC, nZ);
        correctionFactor = tempMean ./ depthMean;
end

% Apply correction factors (except histogram match which already did it)
% to every T-frame uniformly, per Z-plane.
if ~strcmp(p.method, 'histogram match')
    for iC = 1:nC
        if channels(iC)
            for iZ = 1:nZ
                temp(:,:,iC,iZ,:) = temp(:,:,iC,iZ,:) .* correctionFactor(iC,iZ);
            end
        end
    end
end

% Cast back to original class
switch origClass
    case 'uint8'
        imOut = im2uint8(temp);
    case 'uint16'
        imOut = im2uint16(temp);
    otherwise
        imOut = im2single(temp);
end

results.correctionFactor = correctionFactor;
results.norm             = norm;
results.fitObjects       = fitObjects;
results.gof              = gof;
results.normRefWarning   = normRefWarning;
end
