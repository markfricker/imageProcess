function imOut = correctBackground(imIn, p, bgStats)
%CORRECTBACKGROUND  Remove image background using one of several methods.
%
%   imOut = correctBackground(imIn, p, bgStats)
%
% Applies background correction to each selected channel/Z/T slice.
% Negative values after subtraction are clamped to zero.
%
% METHODS
%   'none'       – pass-through
%   'subtract'   – subtract scalar mean from bgStats (single/each T/mean T)
%   'opening'    – morphological top-hat (imtophat)
%   'reconstruct'– top-hat by reconstruction: image minus its opening by
%                  reconstruction (imerode + imreconstruct)
%   'sub low pass'– open → low-pass → subtract (pad/unpad)
%   'surface fit'– fit surface to regional minima, subtract (requires
%                  Curve Fitting Toolbox)
%
% INPUTS
%   imIn    – [nY nX nC nZ nT] numeric array
%   p       – parameter struct:
%               p.method         – method string (see above)
%               p.filterSize     – disk radius for morphological methods
%               p.surfaceModel   – fit model string for 'surface fit'
%                                  (e.g. 'linearinterp','cubicinterp')
%               p.subtractMethod – 'single'|'each T'|'mean T'
%               p.channels       – [1×nC] logical; which channels to correct
%   bgStats – struct from measureBackground with fields .mean [nC×nT]
%             (only needed for 'subtract'; pass [] for other methods)
%
% OUTPUT
%   imOut   – corrected array, same size and class as imIn

[~, ~, nC, nZ, nT] = size(imIn);
imOut    = imIn;
channels = logical(p.channels(1:nC));

if strcmp(p.method, 'none')
    return
end

r  = double(p.filterSize);
SE = strel('disk', r);

for iT = 1:nT
    for iZ = 1:nZ
        for iC = 1:nC
            if ~channels(iC)
                continue
            end

            im   = imIn(:,:,iC,iZ,iT);
            temp = zeros(size(im), 'like', im);

            switch p.method
                case 'subtract'
                    if ~isempty(bgStats) && isfield(bgStats,'mean') && ~isempty(bgStats.mean)
                        switch p.subtractMethod
                            case 'single'
                                mn = bgStats.mean(iC, 1);
                            case 'each T'
                                mn = bgStats.mean(iC, iT);
                            case 'mean T'
                                mn = mean(bgStats.mean(iC,:));
                        end
                        if isinteger(im)
                            temp = imsubtract(im, mn);
                        else
                            temp = im - mn;
                        end
                    else
                        temp = im;   % no ROI measured — pass through
                    end

                case 'opening'
                    temp = imtophat(im, SE);

                case 'reconstruct'
                    % top-hat by reconstruction: the opening by
                    % reconstruction (erode, then grow back under the
                    % image) keeps the background and every structure
                    % wider than the disk, with its true shape; subtracting
                    % it leaves the narrower structures. (Before 2026-09-28
                    % this returned the reconstruction itself -- the
                    % background estimate, not the corrected image -- as
                    % AnalyzER v2 did.)
                    Io   = imerode(im, SE);
                    bg   = imreconstruct(Io, im);
                    if isinteger(im)
                        temp = imsubtract(im, bg);
                    else
                        temp = im - bg;
                    end

                case 'sub low pass'
                    im2  = padarray(im, [r r], 'replicate');
                    im3  = imopen(im2, SE);
                    temp = im2 - imgaussfilt(im3, r, 'padding', 'replicate');
                    temp = temp(r+1:end-r, r+1:end-r);

                case 'surface fit'
                    im2    = padarray(im, [r r], 'replicate');
                    switch p.surfaceModel
                        case {'linearinterp','cubicinterp','biharmonicinterp'}
                            npoints = 10000;
                        otherwise
                            npoints = 1000;
                    end
                    regionalmin = bwulterode(imregionalmin(im2));
                    [yp, xp, zp] = find(double(im2) .* regionalmin);
                    qt      = quantile(zp, [0.1 0.9]);
                    exclude = zp < qt(1) | zp > qt(2);
                    inc     = max(1, round(numel(zp) / npoints));
                    fs      = fit([xp(1:inc:end), yp(1:inc:end)], ...
                                  double(zp(1:inc:end)), p.surfaceModel, ...
                                  'exclude', exclude(1:inc:end));
                    [xg, yg] = meshgrid(1:size(im2,2), 1:size(im2,1));
                    back     = fs(xg, yg);
                    temp     = double(im2) - back;
                    temp     = cast(temp(r+1:end-r, r+1:end-r), 'like', im);
            end

            temp(temp < 0) = 0;
            imOut(:,:,iC,iZ,iT) = temp;
        end
    end
end
end
