function tforms = estimateRegistration(imIn, p)
%ESTIMATEREGISTRATION  Estimate per-frame geometric transforms for time-series registration.
%
%   tforms = estimateRegistration(imIn, p)
%
% Computes a geometric transform for each frame relative to a fixed
% reference frame using SURF feature matching, MSER feature matching,
% or intensity-based (imregtform) optimisation.
%
% INPUTS
%   imIn – [nY nX nC nZ nT] numeric; input image stack.
%   p    – parameter struct:
%            p.roi            – [1×4] [x y w h]; ROI rectangle used for
%                               both fixed and moving crops (pixels).
%            p.channel        – integer; channel index used for estimation.
%            p.transform      – 'translation' | 'rigid' | 'similarity' |
%                               'affine' | 'SURF' | 'MSER'
%                               'SURF'/'MSER' use feature matching;
%                               others use imregtform (intensity-based).
%            p.reference      – 'first' | 'previous'; which frame is the
%                               fixed reference for each moving frame.
%                               With 'previous' the frame-to-frame
%                               transforms are chained, so the returned
%                               transform still maps each frame onto
%                               frame 1 (use it when the sample changes
%                               too much for every frame to match frame 1).
%                               Any other transform name (e.g. 'manual')
%                               is an error.
%            p.smoothUse      – logical; Gaussian pre-filter before matching.
%            p.smooth         – scalar; sigma for imgaussfilt pre-filter.
%            p.gradientUse    – logical; use image gradient magnitude instead
%                               of raw intensity.
%            p.constraintsUse – logical; tighter optimiser constraints (fewer
%                               iterations) for intensity-based methods.
%            p.phaseEstimateUse – logical; seed intensity optimisation with
%                               a phase-correlation translation estimate.
%
% OUTPUT
%   tforms – {1 × nT} cell array; each cell holds a geometric transform
%            object (affine2d or similar) mapping that frame to the fixed
%            reference.  Frame 1 is always the identity transform.
%
% NOTE ON TRANSFORM NAMING
%   'SURF' and 'MSER' use estimateGeometricTransform (similarity transform).
%   'translation', 'rigid', 'similarity', 'affine' use imregtform.
%   Set p.transform = 'rigid' for standard time-lapse drift correction.

[nY, nX, ~, ~, nT] = size(imIn); %#ok<ASGLU>

rect    = p.roi;
channel = p.channel;
method  = p.transform;

supported = {'translation', 'rigid', 'similarity', 'affine', 'SURF', 'MSER'};
if ~any(strcmp(method, supported))
    error('estimateRegistration:unsupportedTransform', ...
        'Transform ''%s'' is not supported; use one of: %s.', method, strjoin(supported, ', '));
end
previousRef = strcmpi(p.reference, 'previous');

tforms      = cell(1, nT);
tforms{1}   = affine2d;          % identity for first frame

% Pre-compute fixed reference image (frame 1, cropped)
fixedRaw = squeeze(imIn( ...
    rect(2) : rect(2)+rect(4), ...
    rect(1) : rect(1)+rect(3), ...
    channel, 1, 1));

fixedRefObj = imref2d( ...
    [rect(4)+1, rect(3)+1], ...
    [rect(2), rect(2)+rect(4)+1], ...
    [rect(1), rect(1)+rect(3)+1]);

for iT = 2:nT

    % Update fixed reference when using 'previous' mode
    if previousRef && iT > 2
        fixedRaw = squeeze(imIn( ...
            rect(2) : rect(2)+rect(4), ...
            rect(1) : rect(1)+rect(3), ...
            channel, 1, iT-1));
    end

    fixed  = prepareFrame(fixedRaw,  p);
    moving = prepareFrame(squeeze(imIn( ...
        rect(2) : rect(2)+rect(4), ...
        rect(1) : rect(1)+rect(3), ...
        channel, 1, iT)), p);

    switch method

        case 'SURF'
            fixedPts  = detectSURFFeatures(fixed,  'MetricThreshold', 750, 'NumOctaves', 3, 'NumScaleLevels', 5);
            movingPts = detectSURFFeatures(moving, 'MetricThreshold', 750, 'NumOctaves', 3, 'NumScaleLevels', 5);
            [fixedFeat,  fixedValid]  = extractFeatures(fixed,  fixedPts,  'Upright', false);
            [movingFeat, movingValid] = extractFeatures(moving, movingPts, 'Upright', false);
            pairs = matchFeatures(fixedFeat, movingFeat, 'MatchThreshold', 25, 'MaxRatio', 0.25);
            try
                tform = estimateGeometricTransform( ...
                    movingValid(pairs(:,2)), fixedValid(pairs(:,1)), 'similarity');
            catch
                tform = affine2d;
            end

        case 'MSER'
            fixedPts  = detectMSERFeatures(fixed,  'ThresholdDelta', 2.4, 'RegionAreaRange', [20 21000], 'MaxAreaVariation', 0.55);
            movingPts = detectMSERFeatures(moving, 'ThresholdDelta', 2.4, 'RegionAreaRange', [20 21000], 'MaxAreaVariation', 0.55);
            [fixedFeat,  fixedValid]  = extractFeatures(fixed,  fixedPts,  'Upright', false);
            [movingFeat, movingValid] = extractFeatures(moving, movingPts, 'Upright', false);
            pairs = matchFeatures(fixedFeat, movingFeat, 'MatchThreshold', 50, 'MaxRatio', 0.5);
            try
                tform = estimateGeometricTransform( ...
                    movingValid(pairs(:,2)), fixedValid(pairs(:,1)), 'similarity');
            catch
                tform = affine2d;
            end

        otherwise   % intensity-based: translation / rigid / similarity / affine
            [optimizer, metric] = imregconfig('monomodal');
            if p.constraintsUse
                optimizer.MaximumIterations = 100;
                optimizer.RelaxationFactor  = 0.5;
                optimizer.MaximumStepLength = 0.0625;
            else
                optimizer.MaximumIterations = 300;
                optimizer.RelaxationFactor  = 0.75;
                optimizer.MaximumStepLength = 0.1;
            end
            if p.phaseEstimateUse
                tformInit = imregcorr(moving, fixed, 'translation');
                tform = imregtform(moving, fixedRefObj, fixed, fixedRefObj, method, ...
                    optimizer, metric, 'PyramidLevels', 3, 'InitialTransformation', tformInit);
            else
                tform = imregtform(moving, fixedRefObj, fixed, fixedRefObj, method, ...
                    optimizer, metric, 'PyramidLevels', 3);
            end
    end

    % 'previous': tform maps frame iT onto frame iT-1 only. Chain it with the
    % transform of frame iT-1 (which maps that frame onto frame 1) so every
    % returned transform maps its frame onto frame 1, as registerImage
    % applies each one independently. Without this the frame-to-frame
    % shifts were never accumulated and slow drift was not removed.
    if previousRef && iT > 2
        tform = affinetform2d(tformMatrix(tforms{iT-1}) * tformMatrix(tform));
    end

    tforms{iT} = tform;
end

end % estimateRegistration


% =========================================================================
function A = tformMatrix(tf)
%TFORMMATRIX  3x3 matrix of a 2-D geometric transform, premultiply convention
% (x' = A*[x; y; 1]), for both the current objects (affinetform2d,
% rigidtform2d, simtform2d, transltform2d: property A) and legacy affine2d
% (property T, postmultiply).
if isprop(tf, 'A')
    A = tf.A;
else
    A = tf.T';
end
end


% =========================================================================
function imOut = prepareFrame(imIn, p)
%PREPAREFRAME  Optional Gaussian smoothing and gradient magnitude pre-filter.
imOut = imadjust(imIn);
if p.smoothUse
    imOut = imgaussfilt(imOut, p.smooth);
end
if p.gradientUse
    imOut = imgradient(imOut);
end
end
