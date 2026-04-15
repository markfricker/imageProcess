classdef testResampleImage < matlab.unittest.TestCase
%TESTRESAMPLEIMAGE  Unit tests for resampleImage.
%
% Run with:
%   results = runtests('testResampleImage');
%   disp(results)

    methods (Test)

        function testFactor2OutputSize(tc)
            % 64x64x2x3x2 single, factor=2 -> [32 32 2 3 2]
            rng(1);
            imIn = rand(64, 64, 2, 3, 2, 'single');
            p.factor = 2;
            imOut = resampleImage(imIn, p);
            tc.verifySize(imOut, [32 32 2 3 2], ...
                'factor=2 should halve spatial dimensions.');
        end

        function testFactor4OutputSize(tc)
            % 64x64x1x1x1 single, factor=4 -> [16 16 1 1 1]
            rng(2);
            imIn = rand(64, 64, 1, 1, 1, 'single');
            p.factor = 4;
            imOut = resampleImage(imIn, p);
            tc.verifySize(imOut, [16 16 1 1 1], ...
                'factor=4 should reduce 64 to 16 in each spatial dimension.');
        end

        function testLogicalNearestNeighbour(tc)
            % Logical input -> logical output, correct size
            rng(3);
            imIn = rand(64, 64, 1, 1, 1) > 0.5;   % logical
            p.factor = 2;
            imOut = resampleImage(imIn, p);
            tc.verifyClass(imOut, 'logical', ...
                'Logical input should produce logical output.');
            tc.verifySize(imOut, [32 32 1 1 1], ...
                'Logical input factor=2 should give [32 32 1 1 1].');
        end

        function testOutputClassPreserved(tc)
            % uint16 input -> uint16 output
            rng(4);
            imIn = uint16(randi([0 65535], 64, 64, 1, 1, 1));
            p.factor = 2;
            imOut = resampleImage(imIn, p);
            tc.verifyClass(imOut, 'uint16', ...
                'Output class should match input class (uint16).');
        end

        function testFactor1Passthrough(tc)
            % factor=1 -> output same size as input (values may differ
            % slightly due to guided-filter pass being skipped when factor<=1)
            rng(5);
            imIn = rand(32, 32, 2, 2, 2, 'single');
            p.factor = 1;
            imOut = resampleImage(imIn, p);
            tc.verifySize(imOut, size(imIn), ...
                'factor=1 should leave spatial size unchanged.');
        end

    end % methods (Test)
end
