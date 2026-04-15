classdef testProjectImage < matlab.unittest.TestCase
%TESTPROJECTIMAGE  Unit tests for projectImage.
%
% Run with:
%   results = runtests('testProjectImage');
%   disp(results)

    properties
        im5D   % [16 16 2 5 3] single test volume
        pBase  % default parameter struct
    end

    methods (TestMethodSetup)
        function buildTestData(tc)
            rng(42);
            tc.im5D = rand(16, 16, 2, 5, 3, 'single');

            p.method         = 'max';
            p.preFilterSize  = 3;
            p.singleChannels = [false false];
            p.singlePlane    = 1;
            p.edfChannels    = [false false];
            p.edfMethod      = 'max';
            p.inputChannels  = [true true];
            p.applyChannels  = [true true];
            tc.pBase = p;
        end
    end

    % ------------------------------------------------------------------
    methods (Test)

        function testOutputHasSingleZDim(tc)
            % Any projection should collapse nZ to 1
            imOut = projectImage(tc.im5D, tc.pBase);
            tc.verifyEqual(size(imOut, 4), 1, ...
                'Z dimension of output should be 1 after projection.');
        end

        function testMaxProjection(tc)
            % method='max' -> compare to max(im5D,[],4)
            p = tc.pBase;
            p.method = 'max';
            imOut = projectImage(tc.im5D, p);
            expected = max(tc.im5D, [], 4);
            tc.verifyEqual(imOut, expected, 'AbsTol', single(1e-5), ...
                'Max projection should match max(im5D,[],4).');
        end

        function testMinProjection(tc)
            % method='min' -> compare to min(im5D,[],4)
            p = tc.pBase;
            p.method = 'min';
            imOut = projectImage(tc.im5D, p);
            expected = min(tc.im5D, [], 4);
            tc.verifyEqual(imOut, expected, 'AbsTol', single(1e-5), ...
                'Min projection should match min(im5D,[],4).');
        end

        function testMeanProjection(tc)
            % method='mean' -> compare to mean(im5D,4)
            p = tc.pBase;
            p.method = 'mean';
            imOut = projectImage(tc.im5D, p);
            expected = mean(tc.im5D, 4);
            tc.verifyEqual(imOut, expected, 'AbsTol', single(1e-5), ...
                'Mean projection should match mean(im5D,4).');
        end

        function testSinglePlaneChannel(tc)
            % singleChannels=[true false], singlePlane=3, method='none'
            % -> channel 1 of output equals Z-plane 3 of input channel 1
            p = tc.pBase;
            p.method         = 'none';
            p.singleChannels = [true false];
            p.singlePlane    = 3;
            imOut = projectImage(tc.im5D, p);
            expected = tc.im5D(:,:,1,3,:);  % [16 16 1 1 3]
            tc.verifyEqual(imOut(:,:,1,1,:), expected, 'AbsTol', single(1e-5), ...
                'Single-plane channel 1 should equal Z=3 of input channel 1.');
        end

        function testOutputClassPreserved(tc)
            % Output class should match input class
            imOut = projectImage(tc.im5D, tc.pBase);
            tc.verifyClass(imOut, class(tc.im5D), ...
                'Output class should match input class.');
        end

    end % methods (Test)
end
