classdef testAlignChannels < matlab.unittest.TestCase
%TESTALIGNCHANNELS  Unit tests for alignChannels and estimateChannelOffsets.
%
% Run with:
%   results = runtests('testAlignChannels');
%   disp(results)

    methods (Test)

        function testZeroOffsetsNoChange(tc)
            % Zero offsets, all channels masked -> output equals input
            rng(20);
            imIn = rand(32, 32, 3, 2, 2, 'single');
            offsets.horizontal = [0 0 0];
            offsets.vertical   = [0 0 0];
            offsets.channels   = logical([1 1 1]);
            imOut = alignChannels(imIn, offsets);
            tc.verifyEqual(imOut, imIn, 'AbsTol', single(1e-5), ...
                'Zero offsets should leave image unchanged.');
        end

        function testOutputSizePreserved(tc)
            % Output size must match input size
            rng(21);
            imIn = rand(32, 32, 3, 2, 2, 'single');
            offsets.horizontal = [5 -3 2];
            offsets.vertical   = [-2 4 -1];
            offsets.channels   = logical([1 1 1]);
            imOut = alignChannels(imIn, offsets);
            tc.verifySize(imOut, size(imIn), ...
                'alignChannels should preserve array size.');
        end

        function testUnmaskedChannelUnchanged(tc)
            % channels=[0 1 1], large offsets -> channel 1 (unmasked) unchanged
            rng(22);
            imIn = rand(32, 32, 3, 1, 1, 'single');
            offsets.horizontal = [0  20  20];
            offsets.vertical   = [0 -20 -20];
            offsets.channels   = logical([0 1 1]);
            imOut = alignChannels(imIn, offsets);
            % Channel 1 must be pixel-identical (no warp applied)
            tc.verifyEqual(imOut(:,:,1,:,:), imIn(:,:,1,:,:), ...
                'AbsTol', single(1e-5), ...
                'Unmasked channel should be pixel-identical to input.');
        end

        function testEstimateRecoversTrueOffset(tc)
            % ch2 is ch1 displaced +8 px right, -5 px up (ty=-5 → upward shift).
            % imregtform returns the CORRECTION transform (moving→fixed), which
            % is the opposite sign: horizontal≈-8, vertical≈+5.
            % alignChannels then applies this correction to realign ch2 onto ch1.
            sz = [64, 64];
            ch1s = single(generateTubuleNetwork(sz, 20, 1.5, 0.02, 99));

            % Displace ch2 right 8, up 5 (affine2d tx=+8, ty=-5)
            T = eye(3);
            T(3,1) =  8;
            T(3,2) = -5;
            ch2s = imwarp(ch1s, affine2d(T), 'OutputView', imref2d(sz), ...
                          'SmoothEdges', true);

            imIn = zeros(64, 64, 2, 1, 1, 'single');
            imIn(:,:,1,1,1) = ch1s;
            imIn(:,:,2,1,1) = ch2s;

            offsets = estimateChannelOffsets(imIn, 1, logical([0 1]));

            % Correction is the inverse of the displacement: h≈-8, v≈+5
            tc.verifyLessThanOrEqual(abs(offsets.horizontal(2) - (-8)), 2, ...
                'Estimated horizontal correction for ch2 should be within 2 px of -8.');
            tc.verifyLessThanOrEqual(abs(offsets.vertical(2) - 5), 2, ...
                'Estimated vertical correction for ch2 should be within 2 px of +5.');
        end

    end % methods (Test)
end
