function im = generateBlobImage(sz, nBlobs, radius, noiseLevel, seed)
%GENERATEBLOBIMAGE  Synthetic 2-D Gaussian-blob image for testing.
%
%   im = generateBlobImage(sz, nBlobs, radius, noiseLevel)
%   im = generateBlobImage(sz, nBlobs, radius, noiseLevel, seed)
%
% Places nBlobs Gaussian blobs of the given radius at random positions on a
% canvas of size sz, adds uniform noise (noiseLevel), and normalises the
% result to [0, 1].
%
% INPUTS
%   sz         – [nY nX] canvas size
%   nBlobs     – number of Gaussian blobs to place
%   radius     – standard deviation (pixels) of each Gaussian blob
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

% Pre-compute coordinate grids
[XX, YY] = meshgrid(1:nX, 1:nY);

% Place each blob
for k = 1:nBlobs
    cx = 1 + rand() * (nX - 1);
    cy = 1 + rand() * (nY - 1);
    blob = exp(-((XX - cx).^2 + (YY - cy).^2) / (2 * radius^2));
    canvas = canvas + blob;
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
