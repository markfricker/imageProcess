classdef testCorrectDepth < matlab.unittest.TestCase
%TESTCORRECTDEPTH  Unit tests for estimateDepth and correctDepth.
%
% Run with:
%   results = runtests('testCorrectDepth');
%   disp(results)

    methods (Test)

        function testCorrectionFlattensWithoutSaturating(tc)
            % An object dimming exponentially with depth comes out flat,
            % and the deep planes are not driven into saturation (the
            % scaling must keep the planes' relative brightness).
            [im, obj] = attenuatedStack(10, 6);
            [dm, dw] = estimateDepth(im, [], 'foreground');
            [out, res] = correctDepth(im, dm, ...
                struct('method', 'exponential', 'norm2Start', true, 'channels', true), dw);
            objMean = planeMeans(out, obj);
            tc.verifyLessThan(max(abs(objMean / objMean(1) - 1)), 0.05, ...
                'The object should be equally bright at every depth after correction.');
            pl = double(out(:, :, 1, 6, 1)) / 65535;
            % (scaled-up noise clips a few pixels; scaling each plane on
            % its own clipped the whole object)
            tc.verifyLessThan(mean(pl(obj) >= 0.999), 0.25, ...
                'A mid-depth plane should not saturate after correction.');
            tc.verifyGreaterThan(res.correctionFactor(end), 2, ...
                'The deepest plane should be scaled up.');
        end

        function testRatioMatchesMeanCurve(tc)
            % With 'mean' mode and the ratio method the factors are the
            % inverse of the measured, normalised curve.
            [im, ~] = attenuatedStack(6, 5);
            dm = estimateDepth(im, [], 'mean');
            [~, res] = correctDepth(im, dm, struct('method', 'ratio', 'norm2Start', true, 'channels', true));
            tc.verifyEqual(res.correctionFactor, dm(1) ./ dm, 'RelTol', 1e-10);
        end

        function testPerSectionFollowsDimmingOfTexturedObject(tc)
            % An object with a wide intensity range: the whole-stack
            % threshold keeps only the brightest pixels in deep planes and
            % under-states the dimming; a threshold per section follows it.
            nZ = 10; lambda = 6;
            truth = exp(-((1:nZ) - 1) / lambda);
            im = texturedStack(nZ, truth, 1:nZ);
            fit = struct('method', 'exponential', 'norm2Start', true, 'channels', true);
            [mW, wW] = estimateDepth(im, [], 'foreground');
            [mS, wS] = estimateDepth(im, [], 'foreground per section');
            [~, rW] = correctDepth(im, mW, fit, wW);
            [~, rS] = correctDepth(im, mS, fit, wS);
            errW = abs(1 / rW.correctionFactor(end) - truth(end));
            errS = abs(1 / rS.correctionFactor(end) - truth(end));
            tc.verifyLessThan(errS, 0.05, 'Per-section fit should end near the true dimming.');
            tc.verifyLessThan(errS, errW, 'Per-section should beat the whole-stack threshold here.');
        end

        function testPerSectionIgnoresEmptySections(tc)
            % Sections past the object are noise: near-zero weight, so
            % the fit still follows the object's dimming.
            nZ = 10; lambda = 6;
            truth = exp(-((1:nZ) - 1) / lambda);
            im = texturedStack(nZ, truth, 1:7);
            [mS, wS] = estimateDepth(im, [], 'foreground per section');
            tc.verifyLessThan(max(wS(8:10)), 1, 'Empty sections should get ~zero weight.');
            tc.verifyGreaterThan(min(wS(1:7)), 100);
            [~, rS] = correctDepth(im, mS, struct('method', 'exponential', 'norm2Start', true, 'channels', true), wS);
            tc.verifyLessThan(abs(1 / rS.correctionFactor(end) - truth(end)), 0.05);
        end

        function testSinglePlaneIsNoOp(tc)
            im = uint16(65535 * rand(16, 16, 2, 1, 3));
            dm = estimateDepth(im, []);
            out = correctDepth(im, dm, struct('method', 'exponential', 'norm2Start', true, 'channels', [true true]));
            tc.verifyEqual(out, im);
        end

    end
end

function [im, obj] = attenuatedStack(nZ, lambda)
rng(1);
obj = mat2gray(imgaussfilt(rand(96, 96), 2)) > 0.55;      % same blobby object in every plane
nT = 2;
im = zeros(96, 96, 1, nZ, nT);
for z = 1:nZ
    for t = 1:nT
        im(:, :, 1, z, t) = 0.02 + 0.8 * exp(-(z - 1) / lambda) * obj + 0.005 * randn(96, 96);
    end
end
im = uint16(65535 * min(1, max(0, im)));
end

function im = texturedStack(nZ, truth, objectZ)
% filaments of widely varying brightness on a dark background, the same in
% every section listed in objectZ, scaled by truth(z); other sections empty
rng(2);
tex = imgaussfilt(rand(96, 96), 1.5);
obj = mat2gray(tex) > 0.6;
bright = obj .* (0.2 + 0.8 * mat2gray(imgaussfilt(rand(96, 96), 4)));   % 0.2..1 across the object
im = zeros(96, 96, 1, nZ, 2);
for z = 1:nZ
    for t = 1:2
        s = 0.02 + 0.005 * randn(96, 96);
        if ismember(z, objectZ), s = s + 0.9 * truth(z) * bright; end
        im(:, :, 1, z, t) = s;
    end
end
im = uint16(65535 * min(1, max(0, im)));
end

function m = planeMeans(im, obj)
nZ = size(im, 4);
m = zeros(1, nZ);
for z = 1:nZ
    pl = double(squeeze(im(:, :, 1, z, :)));
    m(z) = mean(pl(repmat(obj, 1, 1, size(pl, 3))));
end
end
