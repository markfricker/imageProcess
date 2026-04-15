classdef testCropImage < matlab.unittest.TestCase
%TESTCROPIMAGE  Unit tests for cropImage.
%
% Run with:
%   results = runtests('testCropImage');
%   disp(results)

    properties
        im5D  % [20 30 3 5 4] single test volume
    end

    methods (TestMethodSetup)
        function buildTestImage(tc)
            rng(42);
            tc.im5D = rand(20, 30, 3, 5, 4, 'single');
        end
    end

    % ------------------------------------------------------------------
    methods (Test)

        function testNoOpFullRange(tc)
            % roiUse=false, full Z/T range -> output should equal input
            p.roiUse = false;
            p.Zmin   = 1;
            p.Zmax   = 5;
            p.Tmin   = 1;
            p.Tmax   = 4;
            imOut = cropImage(tc.im5D, p, []);
            tc.verifyEqual(imOut, tc.im5D, 'AbsTol', single(1e-5), ...
                'No-op crop should return an array identical to the input.');
        end

        function testZTRangeSubset(tc)
            % No ROI, Z=2:4, T=2:3 -> size [nY nX nC 3 2]
            p.roiUse = false;
            p.Zmin   = 2;
            p.Zmax   = 4;
            p.Tmin   = 2;
            p.Tmax   = 3;
            imOut = cropImage(tc.im5D, p, []);
            tc.verifySize(imOut, [20 30 3 3 2], ...
                'Z/T range subset should produce [20 30 3 3 2].');
            expected = tc.im5D(:,:,:, 2:4, 2:3);
            tc.verifyEqual(imOut, expected, 'AbsTol', single(1e-5), ...
                'Z/T range crop values should match direct indexing.');
        end

        function testAxisAlignedROI(tc)
            % rect=[5 3 10 8], full Z/T -> size [8 10 3 5 4]
            p.roiUse = true;
            p.Zmin   = 1;
            p.Zmax   = 5;
            p.Tmin   = 1;
            p.Tmax   = 4;
            rect = [5 3 10 8];   % [x y w h] -> cols 5:14, rows 3:10
            imOut = cropImage(tc.im5D, p, rect);
            tc.verifySize(imOut, [8 10 3 5 4], ...
                'Axis-aligned ROI should produce [8 10 3 5 4].');
            expected = tc.im5D(3:10, 5:14, :, :, :);
            tc.verifyEqual(imOut, expected, 'AbsTol', single(1e-5), ...
                'Axis-aligned ROI values should match direct indexing.');
        end

        function testZmaxClamped(tc)
            % Zmax=999 should be clamped to nZ=5; output size Z-dim = 5
            p.roiUse = false;
            p.Zmin   = 1;
            p.Zmax   = 999;
            p.Tmin   = 1;
            p.Tmax   = 4;
            imOut = cropImage(tc.im5D, p, []);
            tc.verifySize(imOut, [20 30 3 5 4], ...
                'Zmax=999 should clamp to nZ; output Z-dim should be 5.');
        end

        function testOutputClass(tc)
            % Input is single -> output should remain single
            p.roiUse = false;
            p.Zmin   = 1;
            p.Zmax   = 5;
            p.Tmin   = 1;
            p.Tmax   = 4;
            imOut = cropImage(tc.im5D, p, []);
            tc.verifyClass(imOut, 'single', ...
                'Output class should match input class (single).');
        end

    end % methods (Test)
end
