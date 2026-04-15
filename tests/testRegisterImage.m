classdef testRegisterImage < matlab.unittest.TestCase
%TESTREGISTERIMAGE  Unit tests for registerImage.
%
% Run with:
%   results = runtests('testRegisterImage');
%   disp(results)

    methods (Test)

        function testIdentityPreservesImage(tc)
            % Identity affine2d for all frames -> output equal to input
            rng(10);
            imIn = rand(32, 32, 2, 2, 4, 'single');
            nT = 4;
            tforms = repmat({affine2d()}, 1, nT);
            imOut = registerImage(imIn, tforms);
            tc.verifyEqual(imOut, imIn, 'AbsTol', single(1e-5), ...
                'Identity transform should leave image unchanged.');
        end

        function testOutputSizePreserved(tc)
            % Output size must match input size
            rng(11);
            imIn = rand(32, 48, 2, 3, 4, 'single');
            nT = 4;
            tforms = repmat({affine2d()}, 1, nT);
            imOut = registerImage(imIn, tforms);
            tc.verifySize(imOut, [32 48 2 3 4], ...
                'registerImage should preserve array size.');
        end

        function testOutputClassPreserved(tc)
            % Output class should match input class
            rng(12);
            imIn = uint16(randi([0 1000], 16, 16, 1, 1, 2));
            tforms = {affine2d(), affine2d()};
            imOut = registerImage(imIn, tforms);
            tc.verifyClass(imOut, 'uint16', ...
                'Output class should match input class (uint16).');
        end

        function testKnownTranslation(tc)
            % Gaussian blob at col 10, row 16; translate +5 in x -> peak
            % should move to col 15 (within 1 pixel tolerance).
            nY = 32;  nX = 32;
            [XX, YY] = meshgrid(1:nX, 1:nY);
            blobCX = 10;  blobCY = 16;
            slice = single(exp(-((XX - blobCX).^2 + (YY - blobCY).^2) / (2*2^2)));
            imIn = reshape(slice, nY, nX, 1, 1, 1);

            % affine2d: T(3,1) = dx in x direction
            T = eye(3);
            T(3,1) = 5;   % translate +5 columns
            tform = affine2d(T);

            imOut = registerImage(imIn, {tform});

            % Find peak column of output
            out2D = imOut(:,:,1,1,1);
            [~, linearIdx] = max(out2D(:));
            [~, colPeak] = ind2sub([nY, nX], linearIdx);

            tc.verifyLessThanOrEqual(abs(colPeak - 15), 1, ...
                'Peak column after +5 translation should be within 1 px of column 15.');
        end

    end % methods (Test)
end
