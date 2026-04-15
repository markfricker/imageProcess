function im = generateTubuleNetwork(sz, nTubules, sigma, noiseLevel, seed)
%GENERATETUBULENETWORK  Synthetic 2-D tubule-network image for testing.
%
%   im = generateTubuleNetwork(sz, nTubules, sigma, noiseLevel)
%   im = generateTubuleNetwork(sz, nTubules, sigma, noiseLevel, seed)
%
% Draws nTubules random line segments on a blank canvas of size sz, applies
% Gaussian blur (sigma), adds uniform noise (noiseLevel), and normalises
% the result to [0, 1].
%
% INPUTS
%   sz         – [nY nX] canvas size
%   nTubules   – number of random line segments to draw
%   sigma      – Gaussian blur standard deviation (pixels)
%   noiseLevel – amplitude of additive uniform noise (0 = no noise)
%   seed       – (optional) rng seed for reproducibility
%
% OUTPUT
%   im – [sz(1) sz(2)] double image in [0, 1]

if nargin >= 5
    rng(seed);
end

nY = sz(1);
nX = sz(2);
canvas = zeros(nY, nX);

% Draw each tubule as a rasterised line segment
for k = 1:nTubules
    % Random start and end points (fractional, within image bounds)
    x1 = 1 + rand() * (nX - 1);
    y1 = 1 + rand() * (nY - 1);
    x2 = 1 + rand() * (nX - 1);
    y2 = 1 + rand() * (nY - 1);

    % Number of sample points along the segment
    nPts = ceil(max(abs(x2-x1), abs(y2-y1))) * 2 + 2;
    xs = round(linspace(x1, x2, nPts));
    ys = round(linspace(y1, y2, nPts));

    % Clamp to valid range
    xs = max(1, min(nX, xs));
    ys = max(1, min(nY, ys));

    idx = sub2ind([nY, nX], ys, xs);
    canvas(idx) = 1;
end

% Gaussian blur
if sigma > 0
    canvas = imgaussfilt(canvas, sigma);
end

% Add noise
if noiseLevel > 0
    canvas = canvas + noiseLevel * rand(nY, nX);
end

% Normalise to [0, 1]
mn = min(canvas(:));
mx = max(canvas(:));
if mx > mn
    im = (canvas - mn) / (mx - mn);
else
    im = canvas;
end
end
