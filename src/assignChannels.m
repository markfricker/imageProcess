function imOut = assignChannels(imIn, p)
%ASSIGNCHANNELS  Normalise, reorder and convert precision for a 5-D image.
%
%   imOut = assignChannels(imIn, p)
%
% Combines per-channel intensity normalisation (with optional cross-channel
% range linking), channel reordering, and output precision conversion.
% This is the batch-scriptable equivalent of the GUI assignChannel step.
%
% INPUTS
%   imIn  – [nY nX nC nZ nT] numeric array
%   p     – parameter struct:
%             p.channelMap      – [1×6] integer; input channel i maps to
%                                  output channel channelMap(i).
%                                  0 = exclude this input channel.
%             p.channelNorm     – [1×6] logical; true = rescale this
%                                  channel to the full range of the input
%                                  class (integer) or [0,1] (float).
%             p.channelLink     – [1×6] integer; 0 = independent;
%                                  k = share intensity range with channel k
%                                  (both channels use the union of their
%                                  [min, max] extents before normalising).
%             p.channelIdentity – {1×6} cell of label strings; informational
%                                  only, not used in any computation.
%             p.precision       – '8-bit' | '16-bit' | 'single' | 'double'
%
% OUTPUT
%   imOut – reordered, normalised array in the requested precision class.
%
% NOTES
%   • Channels whose channelMap entry is 0 are excluded from the output.
%   • The output has max(channelMap) channels; gaps (if any) are zero.
%   • Cross-channel linking is symmetric and idempotent: linking channel A
%     to B and B to A both produce the same result.

[nY, nX, nC, nZ, nT] = size(imIn);

cMap  = double(p.channelMap(1:6));
cLink = double(p.channelLink(1:6));
cNorm = logical(p.channelNorm(1:6));

% ---- per-channel intensity range (only channels actually consumed below:
% a channel normalised directly, or linked-to by one that is) -------------
needsRange = false(1, nC);
for iC = 1:min(6, nC)
    if cNorm(iC)
        needsRange(iC) = true;
        iL = cLink(iC);
        if iL >= 1 && iL <= nC
            needsRange(iL) = true;
        end
    end
end
minC = zeros(1, nC, 'double');
maxC = zeros(1, nC, 'double');
for iC = find(needsRange)
    minC(iC) = double(min(imIn(:,:,iC,:,:), [], 'all'));
    maxC(iC) = double(max(imIn(:,:,iC,:,:), [], 'all'));
end

% ---- cross-channel linking (symmetric) ----------------------------------
for iC = 1:min(6, nC)
    iL = cLink(iC);
    if iL >= 1 && iL <= nC && iL ~= iC
        sharedMin = min(minC(iC), minC(iL));
        sharedMax = max(maxC(iC), maxC(iL));
        minC(iC) = sharedMin;  minC(iL) = sharedMin;
        maxC(iC) = sharedMax;  maxC(iL) = sharedMax;
    end
end

% ---- full-range target for this input class -----------------------------
if isinteger(imIn)
    mx = double(intmax(class(imIn)));
else
    mx = 1.0;
end

% ---- per-channel normalisation (preserves input class) ------------------
temp = imIn;
for iC = 1:nC
    if iC <= 6 && cNorm(iC) && (maxC(iC) > minC(iC))
        temp(:,:,iC,:,:) = cast( ...
            rescale(imIn(:,:,iC,:,:), 0, mx, ...
                    'InputMin', minC(iC), 'InputMax', maxC(iC)), ...
            class(imIn));
    end
end

% ---- channel reordering -------------------------------------------------
activeRows = find(cMap(1:min(6, nC)));
if isempty(activeRows)
    error('assignChannels:noChannels', ...
        'No active channels: all channelMap entries are 0.');
end

nOutChan = max(cMap(1:min(6, nC)));
if nOutChan > nC
    error('assignChannels:invalidChannelOrder', ...
        ['channelMap references output channel %d but the image has ' ...
         'only %d channels.'], nOutChan, nC);
end

out = zeros(nY, nX, nOutChan, nZ, nT, 'like', temp);
for k = 1:numel(activeRows)
    iIn  = activeRows(k);
    iOut = cMap(iIn);
    out(:,:,iOut,:,:) = temp(:,:,iIn,:,:);
end

% ---- precision conversion -----------------------------------------------
switch p.precision
    case '8-bit'
        imOut = im2uint8(out);
    case '16-bit'
        imOut = im2uint16(out);
    case 'single'
        imOut = im2single(rescale(out));
    case 'double'
        imOut = rescale(double(out));
    otherwise
        error('assignChannels:unknownPrecision', ...
            ['Unknown precision "%s". ' ...
             'Use ''8-bit'', ''16-bit'', ''single'' or ''double''.'], ...
            p.precision);
end
end
