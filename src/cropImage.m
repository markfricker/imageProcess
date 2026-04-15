function imOut = cropImage(imIn, p, rect)
%CROPIMAGE  Crop a 5-D image array by XY ROI and/or Z/T range.
%
%   imOut = cropImage(imIn, p, rect)
%
% Performs axis-aligned or rotated-rectangle XY cropping, combined with
% Z-slice and T-frame range selection.
%
% INPUTS
%   imIn  – [nY nX nC nZ nT] numeric array
%   p     – parameter struct:
%             p.roiUse  – logical; apply XY rect crop when true
%             p.Zmin    – first Z slice to retain (1-based)
%             p.Zmax    – last  Z slice to retain (1-based, Inf = all)
%             p.Tmin    – first T frame to retain (1-based)
%             p.Tmax    – last  T frame to retain (1-based, Inf = all)
%   rect  – [x y w h] or [x y w h theta]; empty → no XY crop
%
% OUTPUT
%   imOut – cropped array (same class as imIn)

[~, ~, nC, nZ, nT] = size(imIn);
Zmin = max(1,  min(p.Zmin, nZ));
Zmax = min(nZ, max(1, p.Zmax));
Tmin = max(1,  min(p.Tmin, nT));
Tmax = min(nT, max(1, p.Tmax));

if ~isempty(rect) && p.roiUse
    if numel(rect) < 5 || rect(5) == 0 || rect(5) == 360
        % Axis-aligned rectangle
        imOut = imIn(rect(2):rect(2)+rect(4)-1, rect(1):rect(1)+rect(3)-1, ...
                     1:nC, Zmin:Zmax, Tmin:Tmax);
    else
        % Rotated rectangle — warp then Z/T slice
        width  = rect(3);
        height = rect(4);
        theta  = rect(5);
        xcen   = rect(1) + width/2;
        ycen   = rect(2) + height/2;
        rot    = [cosd(theta),  sind(theta);
                 -sind(theta),  cosd(theta)];
        trans      = [width, height]/2 - [xcen, ycen]*rot;
        tform      = rigid2d(rot, trans);
        outputView = imref2d([height, width]);
        imOut = imwarp(imIn(:,:,:,Zmin:Zmax,Tmin:Tmax), tform, ...
                       'Interp', 'linear', 'OutputView', outputView);
    end
else
    imOut = imIn(:,:,:, Zmin:Zmax, Tmin:Tmax);
end
end
