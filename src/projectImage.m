function imOut = projectImage(imIn, p)
%PROJECTIMAGE  Z-project a 5-D image array across the Z dimension.
%
%   imOut = projectImage(imIn, p)
%
% Supports three independent per-channel projection modes, selected by
% logical channel-index vectors in p.  Channels may be assigned to at most
% one mode; unassigned channels use the global method.
%
% INPUTS
%   imIn  – [nY nX nC nZ nT] numeric array
%   p     – parameter struct:
%     p.method         – global method string: 'none' | 'max' | 'min' |
%                        'mean' | 'max plane' | 'max + max plane' |
%                        'mean + max plane'
%     p.preFilterSize  – averaging kernel size used by max-plane methods
%     p.singleChannels – logical [1×nC]; channels to take from one Z plane
%     p.singlePlane    – Z index for single-plane channels
%     p.edfChannels    – logical [1×nC]; channels to project with EDF
%     p.edfMethod      – EDF method: 'tenengrad'|'lapvar'|'sml'|
%                        'wavelet'|'brenner'|'min'|'max'
%     p.inputChannels  – logical [1×nC]; reference channels for max-plane
%     p.applyChannels  – logical [1×nC]; channels that receive max-plane
%                        result (others get their own projection)
%
% OUTPUT
%   imOut – [nY nX nC 1 nT] projected array (same class as imIn)
%
% NOTE
%   EDF methods 'tenengrad','lapvar','sml','wavelet','brenner' require the
%   edof_fuse function (supplied in utils/).

[nY, nX, nC, ~, nT] = size(imIn);

singleC = p.singleChannels(1:nC);
edfC    = p.edfChannels(1:nC);
inputC  = p.inputChannels(1:nC);
applyC  = p.applyChannels(1:nC);

% EDF method 'none' = no EDF: the channels ticked for EDF are projected
% with the global method like the rest (before 2026-09-29 they came out
% blank -- edf 'none' had no case below)
if strcmpi(p.edfMethod, 'none')
    edfC(:) = false;
end
% global method 'none' with no single-plane or EDF channels = no projection
% at all: the stack is returned unchanged (before 2026-09-29 the output was
% a blank single plane)
if strcmpi(p.method, 'none') && ~any(singleC | edfC)
    imOut = imIn;
    return
end
imOut = zeros(nY, nX, nC, 1, nT, 'like', imIn);

% --- Global projection (channels not reserved for single or EDF) ----------
if any(~(singleC | edfC))
    switch p.method
        case 'none'
            % Leave as zeros

        case {'max', 'maximum'}
            imOut(:,:,:,1,:) = max(imIn, [], 4);

        case {'min', 'minimum'}
            imOut(:,:,:,1,:) = min(imIn, [], 4);

        case {'mean', 'average'}
            imOut(:,:,:,1,:) = mean(imIn, 4);

        case {'mx plane', 'max plane', 'mx + mx plan', ...
              'max + max plane', 'mean + max plane'}
            h = fspecial('average', p.preFilterSize);
            for iT = 1:nT
                im  = imIn(:,:,:,:,iT);
                imf = imfilter(im(:,:,inputC,:), h);
                mxc = squeeze(max(imf, [], 3));
                [~, hgt] = max(mxc, [], 3);
                [r, c, z] = find(hgt);
                ch = zeros(size(r));
                for iC = 1:nC
                    ch(:) = iC;
                    idx = sub2ind(size(im), r, c, ch, z);
                    imOut(:,:,iC,1,iT) = reshape(im(idx), nY, nX);
                end
            end
    end

    % Override non-apply channels with their own projection
    if any(~applyC)
        switch p.method
            case 'mean + max plane'
                for iT = 1:nT
                    imOut(:,:,~applyC,1,iT) = mean(imIn(:,:,~applyC,:,iT), 4);
                end
            case {'max + max plane', 'mx + mx plan'}
                for iT = 1:nT
                    imOut(:,:,~applyC,1,iT) = max(imIn(:,:,~applyC,:,iT), [], 4);
                end
        end
    end
end

% --- Single-plane channels ------------------------------------------------
if any(singleC)
    iZ = p.singlePlane;
    for iT = 1:nT
        imOut(:,:,singleC,1,iT) = imIn(:,:,singleC,iZ,iT);
    end
end

% --- EDF channels ---------------------------------------------------------
if any(edfC)
    ref = [];
    for iT = 1:nT
        switch p.edfMethod
            case {'tenengrad', 'lapvar', 'sml', 'wavelet', 'brenner'}
                edofIdx = find(edfC);
                for eC = 1:numel(edofIdx)
                    imOut(:,:,edofIdx(eC),1,iT) = edof_fuse( ...
                        squeeze(imIn(:,:,edofIdx(eC),:,iT)), ...
                        'focus', p.edfMethod);
                end
            case 'min'
                imOut(:,:,edfC,1,iT) = min(imIn(:,:,edfC,:,iT), [], 4);
            case 'max'
                imOut(:,:,edfC,1,iT) = max(imIn(:,:,edfC,:,iT), [], 4);
        end
        % Histogram-match subsequent frames to frame 1
        if iT == 1
            ref = imOut(:,:,edfC,1,1);
        else
            imOut(:,:,edfC,1,iT) = imhistmatch(imOut(:,:,edfC,1,iT), ref);
        end
    end
end
end
