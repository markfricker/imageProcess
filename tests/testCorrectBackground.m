classdef testCorrectBackground < matlab.unittest.TestCase
%TESTCORRECTBACKGROUND  Unit tests for correctBackground.
%
% Run with:
%   results = runtests('testCorrectBackground');
%   disp(results)
%
% Test image: a smooth bright ramp (background) + a thin bright line
% (structure narrower than the disk) + a wide bright square (structure
% wider than the disk).

    methods (Static)
        function [im, line, square] = testImage(cls)
            [xx, ~] = meshgrid(1:96, 1:96);
            bg = 20 + 0.2 * xx;                      % ramp background
            line = false(96); line(40:41, 10:86) = true;
            square = false(96); square(60:85, 60:85) = true;
            im = bg + 100 * line + 60 * square;
            if strcmp(cls, 'uint16')
                im = uint16(im);
            else
                im = single(im / 255);
            end
        end
        function p = params(method)
            p = struct('method', method, 'filterSize', 5, 'surfaceModel', 'poly11', ...
                'subtractMethod', 'single', 'channels', true);
        end
    end

    methods (Test)

        function testReconstructRemovesBackground(tc)
            % the corrected image keeps the thin line and removes the ramp
            for cls = {'uint16', 'single'}
                [im, line, square] = testCorrectBackground.testImage(cls{1});
                out = double(correctBackground(im, testCorrectBackground.params('reconstruct'), []));
                bgPix = ~imdilate(line | square, ones(9));   % the ramp, away from both structures
                tc.verifyLessThan(mean(out(bgPix)), 0.05 * mean(out(line)), ...
                    sprintf('%s: background should be ~0 after reconstruct', cls{1}));
                tc.verifyGreaterThan(mean(out(line)), 0.8 * 100 * scale(cls{1}), ...
                    sprintf('%s: the thin line should survive', cls{1}));
            end
        end

        function testReconstructIsImageMinusReconstruction(tc)
            [im, ~, ~] = testCorrectBackground.testImage('single');
            out = correctBackground(im, testCorrectBackground.params('reconstruct'), []);
            bg  = imreconstruct(imerode(im, strel('disk', 5)), im);
            tc.verifyEqual(out, max(im - bg, 0), 'AbsTol', 1e-6);
        end

        function testOpeningMatchesTophat(tc)
            [im, ~, ~] = testCorrectBackground.testImage('uint16');
            out = correctBackground(im, testCorrectBackground.params('opening'), []);
            tc.verifyEqual(out, imtophat(im, strel('disk', 5)));
        end

        function testNoneIsPassThrough(tc)
            [im, ~, ~] = testCorrectBackground.testImage('uint16');
            tc.verifyEqual(correctBackground(im, testCorrectBackground.params('none'), []), im);
        end

        function testUnselectedChannelUntouched(tc)
            [im, ~, ~] = testCorrectBackground.testImage('uint16');
            im2 = cat(3, im, im);
            p = testCorrectBackground.params('reconstruct');
            p.channels = [true false];
            out = correctBackground(im2, p, []);
            tc.verifyEqual(out(:,:,2), im);
            tc.verifyNotEqual(out(:,:,1), im);
        end

    end
end

function s = scale(cls)
if strcmp(cls, 'uint16'), s = 1; else, s = 1 / 255; end
end
