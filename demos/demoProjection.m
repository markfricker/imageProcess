%DEMOPROJECTION  Demonstrate Z-projection methods on a synthetic 5-Z stack.
%
% Creates a 128x128, 2-channel, 5-Z-plane stack (1 T) in which:
%   Channel 1 (ER tubules) is sharpest at Z=3 (Gaussian blur increases
%             with distance from Z=3).
%   Channel 2 (ER bodies)  is sharpest at Z=2.
%
% Three projections are shown:
%   Row 1 – Max projection
%   Row 2 – Mean projection
%   Row 3 – Single plane Z=3 (same plane for both channels)
%
% NOTE: The projectImage API exposes one singlePlane value shared across
% all channels flagged by singleChannels.  To pick different planes for
% different channels, call projectImage separately per channel or use the
% EDF route.  This demo uses Z=3 for both channels for simplicity.

%% ---- Build synthetic 5-Z stack ----------------------------------------
sz  = [128, 128];
nZ  = 5;
ch1 = zeros([sz, 1, nZ, 1], 'single');
ch2 = zeros([sz, 1, nZ, 1], 'single');

for z = 1:nZ
    % ER tubules: in-focus at Z=3; blur increases as (z-3)^2
    sigma1     = 1 + 0.3 * (z - 3)^2;
    tube       = generateTubuleNetwork(sz, 30, sigma1, 0.03, z);
    ch1(:,:,1,z,1) = single(tube);

    % ER bodies: in-focus at Z=2; blur increases as |z-2|
    sigma2     = 3 + 0.4 * abs(z - 2);
    blob       = generateBlobImage(sz, 12, sigma2, 0.02, z + 10);
    ch2(:,:,1,z,1) = single(blob);
end

% Combine into [128 128 2 5 1]
im5D = cat(3, ch1, ch2);

%% ---- Build projection parameter structs --------------------------------
pBase.preFilterSize  = 3;
pBase.singleChannels = [false false];
pBase.singlePlane    = 3;
pBase.edfChannels    = [false false];
pBase.edfMethod      = 'max';
pBase.inputChannels  = [true true];
pBase.applyChannels  = [true true];

% 1. Max projection
pMax          = pBase;
pMax.method   = 'max';
imMax         = projectImage(im5D, pMax);

% 2. Mean projection
pMean         = pBase;
pMean.method  = 'mean';
imMean        = projectImage(im5D, pMean);

% 3. Single plane Z=3 (both channels)
pSingle                = pBase;
pSingle.method         = 'none';
pSingle.singleChannels = [true true];
pSingle.singlePlane    = 3;
imSingle               = projectImage(im5D, pSingle);

%% ---- Figure: 3 rows x 2 cols ------------------------------------------
methodLabels = {'Max projection', 'Mean projection', 'Single plane (Z=3)'};
projections  = {imMax, imMean, imSingle};
chLabels     = {'Ch 1: ER tubules (in-focus at Z=3)', ...
                'Ch 2: ER bodies (in-focus at Z=2)'};

fig = figure('Name', 'Z-Projection Methods', 'NumberTitle', 'off', ...
             'Units', 'normalized', 'Position', [0.05 0.05 0.80 0.85]);

for iRow = 1:3
    proj = projections{iRow};
    for iC = 1:2
        ax = subplot(3, 2, (iRow-1)*2 + iC, 'Parent', fig);
        imagesc(ax, proj(:,:,iC,1,1));
        colormap(ax, gray);
        axis(ax, 'image', 'off');
        title(ax, sprintf('%s\n%s', methodLabels{iRow}, chLabels{iC}), ...
              'FontSize', 8);
        xlabel(ax, 'X (px)');
        ylabel(ax, 'Y (px)');
    end
end

sgtitle(fig, 'Z-Projection Methods — 128x128 Synthetic Stack (nZ=5)', ...
        'FontSize', 11);
drawnow;
fprintf('demoProjection complete.\n');
fprintf('  NOTE: singlePlane=3 applied to both channels (API limitation).\n');
fprintf('        Ch2 in-focus plane is Z=2; single-plane demo uses Z=3.\n');
