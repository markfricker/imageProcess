function [imOut, results] = correctDepth(imIn, depthMean, p)
%CORRECTDEPTH  Correct depth-dependent intensity attenuation from a Z-stack image.
%
%   [imOut, results] = correctDepth(imIn, depthMean, p)
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
%   imIn       – [nY nX nC nZ nT] numeric array
%   depthMean  – [nC × nZ] double from estimateDepth (one curve per
%                channel, pooled over the whole time series)
%   p          – parameter struct:
%                  p.method     – method string (see above)
%                  p.norm2Start – logical; true = normalise to Z=1
%                                 (shallowest plane), false = normalise to
%                                 the last (deepest) plane
%                  p.channels   – [1×nC] logical; which channels to correct
%
% OUTPUTS
%   imOut   – corrected array, same size as imIn, class determined by
%             original class (uint8 / uint16 / single)
%   results – struct with fields:
%               results.correctionFactor – [nC × nZ] applied per-plane factors
%               results.norm             – [nC × nZ] normalised depth curve
%               results.fitObjects       – {1×nC} cfit objects (exp methods only)
%               results.gof              – {1×nC} goodness-of-fit (exp methods only)
%
% See also: estimateDepth, correctBleach (the time-domain analogue)

[~, ~, nC, nZ, ~] = size(imIn);
channels = logical(p.channels(1:nC));
z        = (1:nZ)';
origClass = class(imIn);

% Single-plane (non-Z-stack) input: depth correction is a true no-op here,
% not just "nothing to correct" -- the fit-based methods need at least a
% couple of Z-planes to fit a decay curve to and would error on one point.
if nZ <= 1
    imOut = imIn;
    results.correctionFactor = ones(nC, 1);
    results.norm             = ones(nC, 1);
    results.fitObjects       = cell(1, nC);
    results.gof              = cell(1, nC);
    return
end

% Normalise depth curve
if p.norm2Start
    normRef = depthMean(:, 1);
else
    normRef = depthMean(:, end);
end
normRef(normRef == 0) = 1;   % guard against zero reference
norm = depthMean ./ normRef;

correctionFactor = ones(nC, nZ);
fitObjects       = cell(1, nC);
gof              = cell(1, nC);

% Working copy as double for mat2gray, every Z-plane. mat2gray's own
% min/max is taken per (channel, Z-plane) across the whole T range,
% pooling the most pixels available while keeping each plane's own
% intensity scale internally consistent across time.
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
            [fitObjects{iC}, gof{iC}] = fit(z, norm(iC,:)', 'exp1');
            correctionFactor(iC,:)    = 1 ./ feval(fitObjects{iC}, z)';
        end

    case 'double exponential'
        for iC = 1:nC
            [fitObjects{iC}, gof{iC}] = fit(z, norm(iC,:)', 'exp2');
            correctionFactor(iC,:)    = 1 ./ feval(fitObjects{iC}, z)';
        end

    case 'histogram match'
        % Each Z-plane is matched to the Z=1 reference (pooled over T),
        % not a single (Z,T) frame.
        ref = squeeze(mean(double(imIn(:,:,:,1,:)), 5));   % [nY nX nC]
        for iZ = 1:nZ
            for iC = 1:nC
                if channels(iC)
                    nT = size(imIn,5);
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
end
