classdef testEstimateRegistration < matlab.unittest.TestCase
    % Tests for estimateRegistration: recovery of a steady drift with both
    % reference modes, and rejection of unsupported transform names.

    properties
        drift = [2 3]   % [x y] px per frame
        nT = 5
    end

    methods (Test)

        function testFirstReferenceRecoversDrift(tc)
            [im, p] = tc.driftingStack();
            p.reference = 'first';
            tforms = estimateRegistration(im, p);
            tc.verifyShifts(tforms);
        end

        function testPreviousReferenceAccumulatesDrift(tc)
            % With 'previous' each frame is matched to the frame before it;
            % the returned transform must still map the frame to frame 1,
            % i.e. the frame-to-frame shifts must be accumulated.
            [im, p] = tc.driftingStack();
            p.reference = 'previous';
            tforms = estimateRegistration(im, p);
            tc.verifyShifts(tforms);
        end

        function testPreviousRegisteredFramesAlign(tc)
            [im, p] = tc.driftingStack();
            p.reference = 'previous';
            reg = registerImage(im, estimateRegistration(im, p));
            ref = double(reg(40:120, 40:120, 1, 1, 1));
            last = double(reg(40:120, 40:120, 1, 1, end));
            tc.verifyGreaterThan(corr(ref(:), last(:)), 0.95);
        end

        function testUnsupportedTransformErrors(tc)
            [im, p] = tc.driftingStack();
            p.transform = 'manual';
            tc.verifyError(@() estimateRegistration(im, p), ...
                'estimateRegistration:unsupportedTransform');
        end
    end

    methods
        function [im, p] = driftingStack(tc)
            rng(1);
            base = imgaussfilt(rand(160), 3);
            im = zeros(160, 160, 1, 1, tc.nT);
            for t = 1:tc.nT
                im(:,:,1,1,t) = imtranslate(base, (t-1) * tc.drift, 'FillValues', mean(base(:)));
            end
            im = im2uint16(mat2gray(im));
            p.roi = [20 20 120 120];
            p.channel = 1;
            p.transform = 'translation';
            p.reference = 'first';
            p.smoothUse = false;
            p.smooth = 1;
            p.gradientUse = false;
            p.constraintsUse = false;
            p.phaseEstimateUse = true;
        end

        function verifyShifts(tc, tforms)
            for t = 2:tc.nT
                A = tformMatrix(tforms{t});
                % the transform maps frame t back onto frame 1: shift = -(t-1)*drift
                tc.verifyEqual(A(1:2, 3)', -(t-1) * tc.drift, 'AbsTol', 0.5, ...
                    sprintf('frame %d', t));
            end
        end
    end
end

function A = tformMatrix(tf)
if isprop(tf, 'A')
    A = tf.A;          % premultiply convention (affinetform2d, transltform2d, ...)
else
    A = tf.T';         % legacy postmultiply (affine2d)
end
end
