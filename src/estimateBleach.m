function bleachMean = estimateBleach(imIn, roi)
%ESTIMATEBLEACH  Measure mean intensity over time to characterise photobleaching.
%
%   bleachMean = estimateBleach(imIn, roi)
%
% Computes the mean intensity within an ROI (or the full frame) for each
% channel and time frame.  Uses Z=1 only (assumes projected or single-Z
% input).  The resulting curve can be passed to correctBleach.
%
% INPUTS
%   imIn       – [nY nX nC nZ nT] numeric array
%   roi        – [x y w h] ROI, or [] to use the full frame
%
% OUTPUT
%   bleachMean – [nC × nT] double; mean intensity per channel per frame
%                Returns [] if the ROI is provided but invalid.

[nY, nX, nC, ~, nT] = size(imIn);

if isempty(roi)
    % Full-frame measurement
    pos = [1, 1, nX, nY];
else
    pos = roi;
end

if isempty(pos) || ~all(pos(1:4))
    bleachMean = [];
    return
end

x = pos(1); y = pos(2); w = pos(3); h = pos(4);
roi3D = double(imIn(y:y+h-1, x:x+w-1, :, 1, :));  % [h w nC 1 nT]

bleachMean = squeeze(mean(roi3D, [1 2 4]));          % [nC × nT]
bleachMean = reshape(bleachMean, nC, nT);
end
