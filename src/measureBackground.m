function bgStats = measureBackground(imIn, rect)
%MEASUREBACKGROUND  Measure background intensity statistics from an ROI.
%
%   bgStats = measureBackground(imIn, rect)
%
% Extracts mean and standard deviation of pixel intensities within the
% specified rectangular ROI for each channel and time frame (Z=1 only,
% assuming a projected or single-Z image).  Returns empty if the ROI is
% invalid.
%
% INPUTS
%   imIn   – [nY nX nC nZ nT] numeric array
%   rect   – [x y w h] ROI (1-based pixel coordinates, axis-aligned)
%            pass empty [] to skip measurement (returns empty bgStats)
%
% OUTPUT
%   bgStats – struct with fields:
%               bgStats.mean – [nC × nT] double; per-channel per-frame mean
%               bgStats.sd   – [nC × nT] double; per-channel per-frame SD
%             Returns [] if rect is empty or any coordinate is zero.

bgStats = [];

if isempty(rect) || ~all(rect(1:4))
    return
end

[~, ~, nC, ~, nT] = size(imIn);
x = rect(1); y = rect(2); w = rect(3); h = rect(4);

roi = double(imIn(y:y+h-1, x:x+w-1, 1:nC, 1, 1:nT));

bgStats.mean = squeeze(mean(roi, [1 2 4]));   % [nC × nT]
bgStats.sd   = squeeze(std(roi, 1, [1 2 4])); % [nC × nT]

% Ensure [nC × nT] shape even when nC=1 or nT=1
bgStats.mean = reshape(bgStats.mean, nC, nT);
bgStats.sd   = reshape(bgStats.sd,   nC, nT);
end
