classdef testAssignChannels < matlab.unittest.TestCase
%TESTASSIGNCHANNELS  Unit tests for assignChannels.
%
%   Run with:
%       results = runtests('testAssignChannels');
%       disp(results)

    % ------------------------------------------------------------------
    % Shared helpers
    % ------------------------------------------------------------------
    methods (Static)
        function p = defaultParams(nC)
            %DEFAULTPARAMS  Identity mapping, no normalisation, 8-bit output.
            p.channelMap      = [1:nC, zeros(1, 6-nC)];
            p.channelNorm     = false(1, 6);
            p.channelLink     = zeros(1, 6);
            p.channelIdentity = repmat({'none'}, 1, 6);
            p.precision       = '8-bit';
        end
    end

    % ------------------------------------------------------------------
    % Identity and reordering
    % ------------------------------------------------------------------
    methods (Test)

        function testIdentityMapping(tc)
            %TESTIDENTITYMAPPING  channelMap=[1,2] passes channels through.
            imIn = uint8(cat(3, ones(4,4)*50, ones(4,4)*100));
            imIn = repmat(imIn, [1,1,1,2,3]);   % nC=2, nZ=2, nT=3
            p    = tc.defaultParams(2);

            imOut = assignChannels(imIn, p);

            tc.verifyEqual(size(imOut), size(imIn));
            tc.verifyEqual(class(imOut), 'uint8');
            tc.verifyEqual(imOut(:,:,1,:,:), imIn(:,:,1,:,:));
            tc.verifyEqual(imOut(:,:,2,:,:), imIn(:,:,2,:,:));
        end

        function testChannelSwap(tc)
            %TESTCHANNELSWAP  channelMap=[2,1] swaps channel 1 and 2.
            imIn = uint8(cat(3, ones(4,4)*10, ones(4,4)*200));
            p    = tc.defaultParams(2);
            p.channelMap = [2, 1, 0, 0, 0, 0];

            imOut = assignChannels(imIn, p);

            tc.verifyEqual(imOut(:,:,1,1,1), imIn(:,:,2,1,1), ...
                'Output ch1 should contain input ch2.');
            tc.verifyEqual(imOut(:,:,2,1,1), imIn(:,:,1,1,1), ...
                'Output ch2 should contain input ch1.');
        end

        function testExcludeChannel(tc)
            %TESTEXCLUDECHANNEL  channelMap(2)=0 drops channel 2.
            imIn = uint8(cat(3, ones(4,4)*50, ones(4,4)*100));
            p    = tc.defaultParams(2);
            p.channelMap = [1, 0, 0, 0, 0, 0];

            imOut = assignChannels(imIn, p);

            tc.verifyEqual(size(imOut, 3), 1, ...
                'Output should have exactly 1 channel.');
            tc.verifyEqual(imOut(:,:,1,1,1), imIn(:,:,1,1,1));
        end

        function testSingleChannel(tc)
            %TESTSINGLECHANNEL  Single-channel image passes through unchanged.
            imIn = uint8(randi([0 255], 8, 8, 1, 1, 5));
            p    = tc.defaultParams(1);

            imOut = assignChannels(imIn, p);

            tc.verifyEqual(imOut, imIn);
        end

    end

    % ------------------------------------------------------------------
    % Normalisation
    % ------------------------------------------------------------------
    methods (Test)

        function testNormalisationStretchesRange(tc)
            %TESTNORMALISATIONSTRETCHES  Low-contrast channel is stretched to [0,255].
            imIn = uint8(cat(3, ones(8,8)*100, ones(8,8)*200));
            imIn(1,1,1,1,1) = 80;   % min of ch1 = 80
            imIn(2,2,1,1,1) = 120;  % max of ch1 = 120
            p = tc.defaultParams(2);
            p.channelNorm(1) = true;

            imOut = assignChannels(imIn, p);

            tc.verifyLessThanOrEqual(double(min(imOut(:,:,1,:,:), [], 'all')), 1, ...
                'Normalised channel min should be near 0.');
            tc.verifyGreaterThanOrEqual(double(max(imOut(:,:,1,:,:), [], 'all')), 254, ...
                'Normalised channel max should be near 255.');
        end

        function testNoNormPreservesValues(tc)
            %TESTNONORMPRESERVESVALUES  channelNorm=false keeps original values.
            imIn = uint8(cat(3, ones(4,4)*50, ones(4,4)*100));
            p    = tc.defaultParams(2);
            % both channels: norm=false (default)

            imOut = assignChannels(imIn, p);

            tc.verifyEqual(imOut(:,:,1,1,1), imIn(:,:,1,1,1));
            tc.verifyEqual(imOut(:,:,2,1,1), imIn(:,:,2,1,1));
        end

        function testNormOnlyAffectsSelectedChannel(tc)
            %TESTNORMSELECTIVECHANNEL  Normalising ch1 must not change ch2.
            imIn = uint8(cat(3, ones(4,4)*100, ones(4,4)*50));
            imIn(1,1,1,1,1) = 80;
            imIn(2,2,1,1,1) = 120;
            p = tc.defaultParams(2);
            p.channelNorm(1) = true;
            % channelNorm(2) stays false

            imOut = assignChannels(imIn, p);

            tc.verifyEqual(imOut(:,:,2,1,1), imIn(:,:,2,1,1), ...
                'Unnormalised channel 2 should be unchanged.');
        end

    end

    % ------------------------------------------------------------------
    % Cross-channel linking
    % ------------------------------------------------------------------
    methods (Test)

        function testLinkSharesRange(tc)
            %TESTLINKSHARESRANGE  Linked channels use the union of their ranges.
            % ch1: [10, 100],  ch2: [50, 200]  → linked range: [10, 200]
            imIn = uint8(10 * ones(4,4,2));
            imIn(1,1,1,1,1) = 100;   % max of ch1
            imIn(1,1,2,1,1) = 200;   % max of ch2
            imIn(2,2,2,1,1) = 50;    % min of ch2

            p = tc.defaultParams(2);
            p.channelNorm(1) = true;
            p.channelNorm(2) = true;
            p.channelLink(1) = 2;    % ch1 linked to ch2

            imOut = assignChannels(imIn, p);

            % Both channels should be normalised to the same [10, 200] range.
            % ch1 max (=100) should map to 100/200 * 255 ≈ 127 (not 255).
            ch1Max = double(max(imOut(:,:,1,:,:), [], 'all'));
            tc.verifyLessThan(ch1Max, 145, ...
                'With linked range [10,200], ch1 max (100) should map to ~127, not 255.');
        end

    end

    % ------------------------------------------------------------------
    % Precision conversion
    % ------------------------------------------------------------------
    methods (Test)

        function testPrecisionUint8(tc)
            imIn  = uint16(cat(3, ones(4,4)*1000, ones(4,4)*2000));
            p     = tc.defaultParams(2);
            p.channelNorm = [true, true, false, false, false, false];
            p.precision   = '8-bit';

            imOut = assignChannels(imIn, p);

            tc.verifyEqual(class(imOut), 'uint8');
        end

        function testPrecisionUint16(tc)
            imIn  = uint8(cat(3, ones(4,4)*50, ones(4,4)*100));
            p     = tc.defaultParams(2);
            p.channelNorm = [true, true, false, false, false, false];
            p.precision   = '16-bit';

            imOut = assignChannels(imIn, p);

            tc.verifyEqual(class(imOut), 'uint16');
        end

        function testPrecisionSingle(tc)
            imIn  = uint8(randi([0 255], 4, 4, 2));
            p     = tc.defaultParams(2);
            p.channelNorm = [true, true, false, false, false, false];
            p.precision   = 'single';

            imOut = assignChannels(imIn, p);

            tc.verifyEqual(class(imOut), 'single');
            tc.verifyGreaterThanOrEqual(min(imOut(:)), single(0));
            tc.verifyLessThanOrEqual(max(imOut(:)), single(1));
        end

        function testPrecisionDouble(tc)
            imIn  = uint8(randi([0 255], 4, 4, 2));
            p     = tc.defaultParams(2);
            p.channelNorm = [true, true, false, false, false, false];
            p.precision   = 'double';

            imOut = assignChannels(imIn, p);

            tc.verifyEqual(class(imOut), 'double');
            tc.verifyGreaterThanOrEqual(min(imOut(:)), 0);
            tc.verifyLessThanOrEqual(max(imOut(:)), 1);
        end

        function testSingleInputSingleOutput(tc)
            %TESTSINGEINPUT  Single-precision float input, single output.
            imIn  = single(rand(4, 4, 2));
            p     = tc.defaultParams(2);
            p.channelNorm = [true, true, false, false, false, false];
            p.precision   = 'single';

            imOut = assignChannels(imIn, p);

            tc.verifyEqual(class(imOut), 'single');
        end

    end

    % ------------------------------------------------------------------
    % Error handling
    % ------------------------------------------------------------------
    methods (Test)

        function testErrorAllZeroMap(tc)
            imIn = uint8(ones(4,4,2));
            p    = tc.defaultParams(2);
            p.channelMap = zeros(1, 6);

            tc.verifyError(@() assignChannels(imIn, p), ...
                'assignChannels:noChannels');
        end

        function testErrorMapExceedsChannels(tc)
            imIn = uint8(ones(4,4,2));
            p    = tc.defaultParams(2);
            p.channelMap = [1, 5, 0, 0, 0, 0];  % output ch 5 > nC=2

            tc.verifyError(@() assignChannels(imIn, p), ...
                'assignChannels:invalidChannelOrder');
        end

        function testErrorUnknownPrecision(tc)
            imIn = uint8(ones(4,4,1));
            p    = tc.defaultParams(1);
            p.precision = 'float32';  % unknown

            tc.verifyError(@() assignChannels(imIn, p), ...
                'assignChannels:unknownPrecision');
        end

    end

end
