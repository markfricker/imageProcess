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

function m = planeMeans(im, obj)
nZ = size(im, 4);
m = zeros(1, nZ);
for z = 1:nZ
    pl = double(squeeze(im(:, :, 1, z, :)));
    m(z) = mean(pl(repmat(obj, 1, 1, size(pl, 3))));
end
end
