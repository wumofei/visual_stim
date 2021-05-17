%{
Linearly offset and scale input video to obtain desired changes in
mean luminance and contrast.

# Required arguments:

`input_video`: char or matrix. If char, interpret as video filename.
Otherwise, interpret as 8-bit grayscale video matrix.

`dMean`: nonnegative scalar. For example, if `dMean = 1.5`, adjust video
to obtain 150% of original mean intensity.

`dContrast`: nonnegative scalar. For example, if `dContrast = 0.75`,
adjust video to obtain 75% of original contrast.

# Output:

If `input_video` is a char or string, output video is written to a new
.avi file in the same directory as `input_video`. If `input_video` is a
matrix, output video is given as a workspace variable.

# Bugs:

Division by zero when desired new contrast equals 1 on lines 69 and 98.
Need to fix algorithm for determining new pixel value limits. Consider
looking through inbuilt MATLAB functions `imadjust.m` and
`stretchlim.m`.

%}

function output_video = video_linstretch(input_video, dMean, dContrast)

% Check inputs.
if nargin < 3
    error('3 inputs required: input video, desired change in mean intensity, and desired change in contrast.')
elseif ~ischar(input_video) && ~isstring(input_video) && ~ismatrix(input_video)
    error('Input video must be given as filename or matrix.')
elseif ~isscalar(dMean) || dMean <= 0
    error('Desired change in mean intensity must be given as a positive scalar.')
elseif ~isscalar(dContrast) || dContrast <= 0
    error('Desired change in contrast must be given as a positive scalar.')
end

if ischar(input_video)
    % Read input video.
    videoObj = VideoReader(input_video);
    % Record min, max, and mean for each frame.
    frameMin = NaN(1,videoObj.NumFrames);
    frameMax = NaN(1,videoObj.NumFrames);
    frameMean = NaN(1,videoObj.NumFrames);
    for i = 1:videoObj.NumFrames
        frame = squeeze(read(videoObj,i));
        frameMin(i) = min(frame(:));
        frameMax(i) = max(frame(:));
        frameMean(i) = mean2(frame);
    end
    % Find min, max, and mean for entire video.
    oldMin = min(frameMin);
    oldMax = max(frameMax);
    oldMean = mean(frameMean);
    if oldMax == oldMin
        error('Input video is uniform.')
    end
    oldContrast = (oldMax - oldMin) / (oldMax + oldMin);
    oldMean_norm = (oldMean - oldMin) / (oldMax - oldMin);
    newMean = oldMean * dMean;
    newContrast = oldContrast * dContrast;
    % Find new range to match desired change in mean intensity and contrast.
    newMin = newMean / oldMean_norm / ((1 + newContrast) / (1 - newContrast) + 1/(oldMean_norm) - 1); % needs fixing. Currently cannot accommodate case where `newContrast = 1`.
    newMax = newMin * (1 + newContrast) / (1 - newContrast);
    % Linearly stretch frames to desired range.
    [~,filename] = fileparts(videoObj.Name);
    mkdir(videoObj.Path, [filename sprintf('_dMean%i_dContrast%i', dMean*100, dContrast*100)]);
    for i = 1:videoObj.NumFrames
        frame = squeeze(read(videoObj,i));
        frame = (frame - oldMin) / (oldMax - oldMin) * (newMax - newMin) + newMin;
        imwrite(frame, fullfile(videoObj.Path, [filename sprintf('_dMean%i_dContrast%i', dMean*100, dContrast*100)], sprintf('frame%i.tif',i)));
    end
    % Combine frames into output video.
    outfile = VideoWriter(fullfile(videoObj.Path, [filename sprintf('_dMean%i_dContrast%i.avi', dMean*100, dContrast*100)]), 'Grayscale AVI');
    outfile.FrameRate = videoObj.FrameRate;
    open(outfile);
    for i = 1:videoObj.NumFrames
        frame = imread(fullfile(videoObj.Path, [filename sprintf('_dMean%i_dContrast%i', dMean*100, dContrast*100)], sprintf('frame%i.tif',i)));
        writeVideo(outfile,frame)
    end
    close(outfile);
else
    % Find original mean and contrast.
    oldMin = min(input_video(:));
    oldMax = max(input_video(:));
    oldContrast = (oldMax - oldMin) / (oldMax + oldMin);
    oldMean = mean2(input_video);
    oldMean_norm = (oldMean - oldMin) / (oldMax - oldMin);
    % Determine new range to match desired change in mean and contrast.
    newMean = oldMean * dMean;
    newContrast = oldContrast * dContrast;
    newMin = newMean / oldMean_norm / ((1 + newContrast) / (1 - newContrast) + 1/(oldMean_norm) - 1); % needs fixing. Currently cannot accommodate case where `newContrast = 1`.
    newMax = newMin * (1 + newContrast) / (1 - newContrast);
    % Linearly stretch video to desired range.
    output_video = (input_video - oldMin) / (oldMax - oldMin) * (newMax - newMin) + newMin;
    output_video = max(0, min(1, output_video));
end