function output_video = video_linstretch(input_filename, dMean, dContrast)

%% Read input video.
input_videoobj = VideoReader(input_filename);
[~,filename,~] = fileparts(input_filename);
input_fps = input_videoobj.Framerate;
clearvars input_filename

% Extract frames from video object.
input_video = im2double(squeeze(read(input_videoobj, [1 input_videoobj.Duration*input_fps])));
clearvars input_videoobj

%% Find original mean and contrast.
oldMin = min(input_video(:));
oldMax = max(input_video(:));
oldContrast = (oldMax - oldMin) / (oldMax + oldMin);

oldMean = mean2(input_video);
oldMean_norm = (oldMean - oldMin) / (oldMax - oldMin);

%% Determine new range to match desired change in mean and contrast.
newMean = oldMean * dMean;
newContrast = oldContrast * dContrast;

newMin = newMean / oldMean_norm / ((1 + newContrast) / (1 - newContrast) + 1/(oldMean_norm) - 1);
newMax = newMin * (1 + newContrast) / (1 - newContrast);

%% Linearly stretch video to desired range.
output_video = (input_video - oldMin) / (oldMax - oldMin) * (newMax - newMin) + newMin;
output_video = max(0, min(1, output_video));

%% Write output video.
output_file = VideoWriter([filename '_linstretch.avi'], 'Grayscale AVI');
output_file.FrameRate = input_fps;
open(output_file);
writeVideo(output_file, output_video);
close(output_file);

end