function [imOut, results] = correctBleach(imIn, bleachMean, p)
%CORRECTBLEACH  Correct photobleaching from a time-series image.
%
%   [imOut, results] = correctBleach(imIn, bleachMean, p)
%
% Estimates per-frame, per-channel correction factors from a bleaching
% curve (bleachMean from estimateBleach) and applies them to the image.
% Output intensities are normalised to [0,1] via mat2gray before scaling.
%
% ONE correction factor per (channel, frame) is applied UNIFORMLY to every
% Z-plane: standard photobleaching is a global photophysical depletion
% over TIME, not a function of depth, so a single curve per channel is the
% right model (a separate depth-dependent intensity correction, if added,
% would have its own per-Z profile instead). Before 2026-10-01 this
% function silently processed Z=1 only AND left every other Z-plane as
% all-zero in the output (temp was pre-allocated over the full size but
% only ever written at index 1) -- i.e. real data loss on a genuine
% Z-stack with Project disabled, not just "unsupported".
%
% METHODS
%   'ratio'            – correction = 1 / normalised_mean(t)
%   'exponential'      – fit single exponential (exp1), correct by 1/fit(t)
%   'double exponential'– fit double exponential (exp2), correct by 1/fit(t)
%   'histogram match'  – match each frame's histogram to frame 1 (no fitting)
%
% REQUIREMENTS
%   'exponential' and 'double exponential' require Curve Fitting Toolbox.
%
% INPUTS
%   imIn       – [nY nX nC nZ nT] numeric array
%   bleachMean – [nC × nT] double from estimateBleach (one curve per
%                channel, already pooled over Z per its own zSource)
%   p          – parameter struct:
%                  p.method     – method string (see above)
%                  p.norm2Start – logical; true = normalise to frame 1,
%                                 false = normalise to last frame
%                  p.channels   – [1×nC] logical; which channels to correct
%
% OUTPUTS
%   imOut   – corrected array, same size as imIn, class determined by
%             original class (uint8 / uint16 / single)
%   results – struct with fields:
%               results.correctionFactor – [nC × nT] applied per-frame factors
%               results.norm             – [nC × nT] normalised bleach curve
%               results.fitObjects       – {1×nC} cfit objects (exp methods only)
%               results.gof              – {1×nC} goodness-of-fit (exp methods only)
%               results.normRefWarning   – '' if fine, else a message naming
%                                          which channel(s) have an unstable
%                                          normalisation reference (see below)
%
% NORMALISATION REFERENCE STABILITY (p.norm2Start)
%   'ratio'/'exponential'/'double exponential' all normalise the curve by
%   dividing by the reference frame's value (frame 1 or last) before
%   fitting/inverting. If that reference frame happens to be close to zero
%   (e.g. the sample hadn't been excited/focused yet, or a dropped frame),
%   every OTHER frame's correctionFactor = normRef/value collapses toward
%   zero, crushing the corrected image everywhere except at the reference
%   frame itself -- not a subtle effect, a "daft multiplier". This is
%   detected (reference < 10% of the channel's own peak) and returned in
%   results.normRefWarning (and raised via warning()) rather than silently
%   producing a badly-scaled correction; it does not block the correction.
%
% See also: estimateBleach, correctDepth (the depth-domain analogue)

[~, ~, nC, nZ, nT] = size(imIn);
channels = logical(p.channels(1:nC));
x        = (1:nT)';
origClass = class(imIn);

% Normalise bleach curve
if p.norm2Start
    normRef = bleachMean(:, 1);
    refLabel = 'first';
else
    normRef = bleachMean(:, end);
    refLabel = 'last';
end
normRefRaw = normRef;   % pre-guard, for the stability check below
normRef(normRef == 0) = 1;   % guard against exact-zero reference (division)
norm = bleachMean ./ normRef;

% See header: a reference frame close to zero (not just exactly zero)
% still makes every other frame's correction factor collapse toward zero
% once inverted. Flag it rather than silently returning a badly-scaled
% correction -- it does not block the method the caller asked for.
normRefWarning = '';
if ~strcmp(p.method, 'histogram match')
    peakVal = max(bleachMean, [], 2);
    smallRefFrac = 0.1;   % reference frame is < 10% of this channel's own peak
    isBadRef = channels(:) & (normRefRaw < smallRefFrac * peakVal);
    if any(isBadRef)
        normRefWarning = sprintf(['Bleach correction: channel(s) %s -- the %s-frame reference value is ' ...
            'less than %.0f%% of that channel''s peak intensity. Normalising to it will make the ' ...
            'correction factor collapse toward zero everywhere else, crushing the corrected image. ' ...
            'This usually means that frame has little/no real signal -- this correction is probably ' ...
            'not appropriate for this curve as configured; consider normalising to the other end instead.'], ...
            mat2str(find(isBadRef)'), refLabel, 100*smallRefFrac);
        warning('correctBleach:unstableNormReference', '%s', normRefWarning);
    end
end

correctionFactor = ones(nC, nT);
fitObjects       = cell(1, nC);
gof              = cell(1, nC);

% Working copy as double for mat2gray, every Z-plane (not just Z=1 --
% see header). mat2gray's own min/max is taken per (channel, Z-plane)
% across the whole T range, pooling the most pixels available while still
% keeping each Z-plane's own intensity scale internally consistent across
% time (matches the original per-channel design, extended per-Z).
temp = zeros(size(imIn));
for iZ = 1:nZ
    for iC = 1:nC
        temp(:,:,iC,iZ,:) = mat2gray(imIn(:,:,iC,iZ,:));
    end
end

switch p.method
    case 'ratio'
        for iC = 1:nC
            correctionFactor(iC,:) = 1 ./ norm(iC,:);
        end

    case 'exponential'
        for iC = 1:nC
            [fitObjects{iC}, gof{iC}] = fit(x, norm(iC,:)', 'exp1');
            correctionFactor(iC,:)    = 1 ./ feval(fitObjects{iC}, x)';
        end

    case 'double exponential'
        for iC = 1:nC
            [fitObjects{iC}, gof{iC}] = fit(x, norm(iC,:)', 'exp2');
            correctionFactor(iC,:)    = 1 ./ feval(fitObjects{iC}, x)';
        end

    case 'histogram match'
        % Each Z-plane is matched to ITS OWN first-frame reference, not a
        % single plane's histogram applied across the whole stack.
        for iZ = 1:nZ
            ref = squeeze(imIn(:,:,:,iZ,1));
            for iC = 1:nC
                if channels(iC)
                    for iT = 1:nT
                        temp(:,:,iC,iZ,iT) = mat2gray( ...
                            imhistmatch(imIn(:,:,iC,iZ,iT), ref(:,:,iC)));
                    end
                end
            end
        end
        tempMean = squeeze(mean(temp, [1 2 4]));    % pooled over space+Z -> [nC × nT]
        tempMean = reshape(tempMean, nC, nT);
        correctionFactor = tempMean ./ bleachMean;
end

% Apply correction factors (except histogram match which already did it)
% to every Z-plane uniformly.
if ~strcmp(p.method, 'histogram match')
    for iC = 1:nC
        if channels(iC)
            cf = reshape(correctionFactor(iC,:), 1, 1, 1, 1, nT);
            for iZ = 1:nZ
                temp(:,:,iC,iZ,:) = temp(:,:,iC,iZ,:) .* cf;
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
