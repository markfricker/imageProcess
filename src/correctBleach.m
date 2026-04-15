function [imOut, results] = correctBleach(imIn, bleachMean, p)
%CORRECTBLEACH  Correct photobleaching from a time-series image.
%
%   [imOut, results] = correctBleach(imIn, bleachMean, p)
%
% Estimates per-frame, per-channel correction factors from a bleaching
% curve (bleachMean from estimateBleach) and applies them to the image.
% Output intensities are normalised to [0,1] via mat2gray before scaling.
% Only Z=1 is processed (assumes projected or single-Z input).
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
%   imIn       – [nY nX nC nZ nT] numeric array (nZ=1 expected)
%   bleachMean – [nC × nT] double from estimateBleach
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

[~, ~, nC, ~, nT] = size(imIn);
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

% Working copy as double for mat2gray
temp = zeros(size(imIn));
for iC = 1:nC
    temp(:,:,iC,1,:) = mat2gray(imIn(:,:,iC,1,:));
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
        ref = squeeze(imIn(:,:,:,1,1));
        for iC = 1:nC
            if channels(iC)
                for iT = 1:nT
                    temp(:,:,iC,1,iT) = mat2gray( ...
                        imhistmatch(imIn(:,:,iC,1,iT), ref(:,:,iC)));
                end
            end
        end
        tempMean = squeeze(mean(temp, [1 2 4]));    % [nC × nT]
        tempMean = reshape(tempMean, nC, nT);
        correctionFactor = tempMean ./ bleachMean;
end

% Apply correction factors (except histogram match which already did it)
if ~strcmp(p.method, 'histogram match')
    for iC = 1:nC
        if channels(iC)
            cf = reshape(correctionFactor(iC,:), 1, 1, 1, 1, nT);
            temp(:,:,iC,1,:) = temp(:,:,iC,1,:) .* cf;
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
