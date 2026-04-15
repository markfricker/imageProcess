%DEMOCHANNELALIGNMENT  Demonstrate channel offset estimation and correction.
%
% Generates a 3-channel 128x128 synthetic image with known translational
% misalignments, calls estimateChannelOffsets, then alignChannels, and
% visualises the result before and after correction.
%
% Channel layout:
%   Ch 1 (reference) : ER tubule network
%   Ch 2             : Ch1 shifted +12 x, +8 y  + ER bodies blended in
%   Ch 3             : Ch1 shifted  -9 x, -6 y  + ER bodies blended in
%
% True offsets (relative to Ch1):
%   Ch 2: horizontal = +12, vertical = +8
%   Ch 3: horizontal =  -9, vertical = -6

%% ---- Build synthetic channels ------------------------------------------
sz = [128, 128];
base = generateTubuleNetwork(sz, 35, 1.5, 0.03, 1);     % reference ER network
base = single(base);

% Warp ch2: +12 in x, +8 in y
T2       = eye(3); T2(3,1) = 12; T2(3,2) = 8;
tform2   = affine2d(T2);
ref2d    = imref2d(sz);
ch2raw   = imwarp(base, tform2, 'OutputView', ref2d, 'SmoothEdges', true);

% Warp ch3: -9 in x, -6 in y
T3       = eye(3); T3(3,1) = -9; T3(3,2) = -6;
tform3   = affine2d(T3);
ch3raw   = imwarp(base, tform3, 'OutputView', ref2d, 'SmoothEdges', true);

% Blend in ER-body component for realism
bodies   = single(generateBlobImage(sz, 8, 4, 0.01, 2));
ch2      = ch2raw + 0.3 * bodies;
ch3      = ch3raw + 0.3 * bodies;

% Normalise blended channels to [0,1]
ch2 = (ch2 - min(ch2(:))) / (max(ch2(:)) - min(ch2(:)) + eps);
ch3 = (ch3 - min(ch3(:))) / (max(ch3(:)) - min(ch3(:)) + eps);

% Pack into 5-D array [nY nX nC nZ nT]
im5D = zeros(128, 128, 3, 1, 1, 'single');
im5D(:,:,1,1,1) = base;
im5D(:,:,2,1,1) = ch2;
im5D(:,:,3,1,1) = ch3;

%% ---- Estimate and apply offsets ----------------------------------------
refChannel  = 1;
channelMask = logical([0 1 1]);   % estimate for ch2 and ch3

offsets = estimateChannelOffsets(im5D, refChannel, channelMask);

fprintf('\n--- Channel offset estimation ---\n');
fprintf('True offsets  : Ch2 horiz=%+.1f  vert=%+.1f  |  Ch3 horiz=%+.1f  vert=%+.1f\n', ...
    12, 8, -9, -6);
fprintf('Estimated     : Ch2 horiz=%+.2f  vert=%+.2f  |  Ch3 horiz=%+.2f  vert=%+.2f\n', ...
    offsets.horizontal(2), offsets.vertical(2), ...
    offsets.horizontal(3), offsets.vertical(3));

imAligned = alignChannels(im5D, offsets);

%% ---- Figure 1: individual channels before / after ---------------------
fig1 = figure('Name', 'Channel Alignment — Individual Channels', ...
               'NumberTitle', 'off', 'Units', 'normalized', ...
               'Position', [0.05 0.1 0.85 0.75]);

channelNames = {'Ch 1 (reference)', 'Ch 2 (shifted +12x, +8y)', 'Ch 3 (shifted -9x, -6y)'};
for iC = 1:3
    % Row 1: before
    ax = subplot(2, 3, iC, 'Parent', fig1);
    imagesc(ax, im5D(:,:,iC,1,1));
    colormap(ax, gray);
    axis(ax, 'image', 'off');
    title(ax, ['Before: ' channelNames{iC}], 'FontSize', 9);

    % Row 2: after
    ax2 = subplot(2, 3, iC + 3, 'Parent', fig1);
    imagesc(ax2, imAligned(:,:,iC,1,1));
    colormap(ax2, gray);
    axis(ax2, 'image', 'off');
    title(ax2, ['After: ' channelNames{iC}], 'FontSize', 9);
end
sgtitle(fig1, 'Channel Alignment — Individual Channels', 'FontSize', 11);

%% ---- Figure 2: RGB composite before / after ---------------------------
fig2 = figure('Name', 'Channel Alignment — RGB Composite', ...
               'NumberTitle', 'off', 'Units', 'normalized', ...
               'Position', [0.1 0.1 0.70 0.50]);

% Before: misaligned channels produce coloured fringing
rgbBefore = cat(3, im5D(:,:,1,1,1), im5D(:,:,2,1,1), im5D(:,:,3,1,1));
rgbBefore = rgbBefore / (max(rgbBefore(:)) + eps);

% After: aligned channels overlap to produce white/grey
rgbAfter = cat(3, imAligned(:,:,1,1,1), imAligned(:,:,2,1,1), imAligned(:,:,3,1,1));
rgbAfter  = rgbAfter / (max(rgbAfter(:)) + eps);

subplot(1, 2, 1);
imshow(rgbBefore);
title('Before alignment (coloured fringing)', 'FontSize', 10);

subplot(1, 2, 2);
imshow(rgbAfter);
title('After alignment (grey/white overlap)', 'FontSize', 10);

sgtitle(fig2, 'RGB Composite: Ch1=R, Ch2=G, Ch3=B', 'FontSize', 11);

drawnow;
fprintf('\nDone. Figures 1 and 2 show individual channels and RGB composite.\n');
