%DEMOREGISTRATION  Demonstrate frame-to-frame registration on a drifting time series.
%
% Creates a 128x128, 1-channel, 12-frame synthetic ER-network time series
% with cumulative linear drift (dx=1.8*(iT-1), dy=0.9*(iT-1)) plus random
% ER dynamics added per frame.
%
% Steps:
%   1. Build drifted time series
%   2. Estimate per-frame translation transforms (imregtform, monomodal)
%   3. Apply transforms with registerImage
%   4. Visualise: first frame / last frame / kymographs before and after

%% ---- Parameters --------------------------------------------------------
nY    = 128;
nX    = 128;
nT    = 12;
kyRow = 64;          % row used for kymograph

%% ---- Build drifted time series -----------------------------------------
fprintf('Building synthetic drifted time series (%d frames)...\n', nT);

% Static base ER network
base = single(generateTubuleNetwork([nY, nX], 40, 1.5, 0.02, 1));

im5D = zeros(nY, nX, 1, 1, nT, 'single');
for iT = 1:nT
    % Dynamic component: random tubule appearance / disappearance
    dynamic = single(generateTubuleNetwork([nY, nX], 8, 1.5, 0.01, iT * 100));
    frame   = base + 0.15 * dynamic;
    % Re-normalise to [0,1]
    frame   = (frame - min(frame(:))) / (max(frame(:)) - min(frame(:)) + eps);

    % Cumulative drift
    dx = 1.8 * (iT - 1);
    dy = 0.9 * (iT - 1);
    T       = eye(3);
    T(3,1)  = dx;
    T(3,2)  = dy;
    tDrift  = affine2d(T);
    ref2d   = imref2d([nY, nX]);
    frame   = imwarp(frame, tDrift, 'OutputView', ref2d, 'SmoothEdges', true);

    im5D(:,:,1,1,iT) = frame;
end

%% ---- Estimate per-frame transforms -------------------------------------
fprintf('Estimating registration transforms...\n');
[optimizer, metric] = imregconfig('monomodal');
fixed  = im5D(:,:,1,1,1);          % frame 1 is the reference

tforms      = cell(1, nT);
tforms{1}   = affine2d();           % identity for frame 1

dx_est = zeros(1, nT);
dy_est = zeros(1, nT);

for iT = 2:nT
    moving       = im5D(:,:,1,1,iT);
    tforms{iT}   = imregtform(moving, fixed, 'translation', optimizer, metric);
    dx_est(iT)   = tforms{iT}.T(3,1);
    dy_est(iT)   = tforms{iT}.T(3,2);
end

fprintf('\n--- Recovered translation offsets per frame ---\n');
fprintf('Frame  True dx   True dy   Est dx   Est dy\n');
for iT = 1:nT
    trueDx = -1.8 * (iT - 1);    % sign: registration corrects the drift
    trueDy = -0.9 * (iT - 1);
    fprintf('  %2d   %+7.1f   %+7.1f   %+7.2f   %+7.2f\n', ...
        iT, trueDx, trueDy, dx_est(iT), dy_est(iT));
end

%% ---- Apply registration ------------------------------------------------
fprintf('\nApplying transforms with registerImage...\n');
imRegistered = registerImage(im5D, tforms);

%% ---- Figure: 2x2 layout ------------------------------------------------
fig = figure('Name', 'Registration Demo', 'NumberTitle', 'off', ...
             'Units', 'normalized', 'Position', [0.05 0.05 0.80 0.75]);

% (1,1) First frame
ax11 = subplot(2, 2, 1, 'Parent', fig);
imagesc(ax11, im5D(:,:,1,1,1));
colormap(ax11, gray);
axis(ax11, 'image');
title(ax11, 'Frame 1 (reference)', 'FontSize', 10);
xlabel(ax11, 'X (px)');  ylabel(ax11, 'Y (px)');

% (1,2) Last frame (unregistered) — shows drift
ax12 = subplot(2, 2, 2, 'Parent', fig);
imagesc(ax12, im5D(:,:,1,1,end));
colormap(ax12, gray);
axis(ax12, 'image');
title(ax12, sprintf('Frame %d — unregistered (drift visible)', nT), 'FontSize', 10);
xlabel(ax12, 'X (px)');  ylabel(ax12, 'Y (px)');

% (2,1) Kymograph before registration
kyBefore = squeeze(im5D(kyRow, :, 1, 1, :));   % [nX, nT]
ax21 = subplot(2, 2, 3, 'Parent', fig);
imagesc(ax21, kyBefore');   % show frames on x-axis, x-position on y-axis
colormap(ax21, gray);
axis(ax21, 'tight');
title(ax21, sprintf('Kymograph (row %d) — before registration', kyRow), 'FontSize', 10);
xlabel(ax21, 'X (px)');  ylabel(ax21, 'Frame');
% Annotate drift arrow
hold(ax21, 'on');
quiver(ax21, 5, nT*0.15, 10, nT*0.65, 0, 'r', 'LineWidth', 2, ...
       'MaxHeadSize', 0.5, 'AutoScale', 'off');
text(ax21, 7, nT*0.85, 'Drift \rightarrow', 'Color', 'r', 'FontSize', 8);
hold(ax21, 'off');

% (2,2) Kymograph after registration
kyAfter = squeeze(imRegistered(kyRow, :, 1, 1, :));
ax22 = subplot(2, 2, 4, 'Parent', fig);
imagesc(ax22, kyAfter');
colormap(ax22, gray);
axis(ax22, 'tight');
title(ax22, sprintf('Kymograph (row %d) — after registration', kyRow), 'FontSize', 10);
xlabel(ax22, 'X (px)');  ylabel(ax22, 'Frame');
hold(ax22, 'on');
text(ax22, nX*0.5, nT*0.1, 'Corrected', 'Color', 'g', ...
     'FontSize', 8, 'HorizontalAlignment', 'center');
hold(ax22, 'off');

sgtitle(fig, 'Frame Registration — Synthetic Drifting ER Network (12 frames)', ...
        'FontSize', 11);
drawnow;
fprintf('\ndemoRegistration complete.\n');
