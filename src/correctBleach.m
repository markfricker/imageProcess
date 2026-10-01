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

[~, ~, nC, nZ, nT] = size(imIn);
channels = logical(p.channels(1:nC));
x        = (1:nT)';
origClass = class(imIn);

% Normalise bleach curve
if p.norm2Start
    normRef = bleachMean(:, 1);
else
    normRef = bleachMean(:, end);
end
normRef(normRef == 0) = 1;   % guard against zero reference
norm = bleachMean ./ normRef;

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
end
