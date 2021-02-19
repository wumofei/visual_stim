% Read video file into memory as matrix.
% Will exceed memory limits for large videos. Designed for Wei lab OLED
% stimulus, i.e. 800x600 pixel 8-bit grayscale videos 10 seconds long @
% 60 fps, equivalent to 288 megabytes.

function [videomat, fps] = readmovie(input_filename)

input_videoobj = VideoReader(input_filename);
fps = input_videoobj.Framerate;
videomat = squeeze(read(input_videoobj));

end