%{
Convert input video to Wei lab rig OLED display stimulus. Match size of
depicted scene to size of target retina and optimal viewing distance.
Convolve with 3x3 kernel to detect edges. Scale video size to match OLED
size: 800x600 pixels. Downsample to match bipolar cell RF diamater to
produce background checkerboard flicker.

Two ways of selecting region of input video to crop are supported:

1.
Specify `input_dist`, `input_focal_length`, and either
`input_sensor_size` or `input_pixel_size`. Along with
`visualAngle2retinaDist`, `output_retina_size` and
`output_viewing_dist`, pixel dimensions of input video that match
requested optimal viewing distance are calculated in lines 208-215. A
pop-up window displaying the first frame of the input video, overlaid
with a rectangular indicator with the calculated crop pixel dimensions,
will appear. Move rectangular indicator to select desired region to
crop. Double-click to finish. You may change the size of the rectangular
indicator, but keep in mind that you might want to keep the aspect ratio
constant, lest the input video be distorted.

2.
Directly specify a rectangular region in pixels to crop using
`cropRect`. Keep in mind that the Wei lab OLED is 800x600 pixels, so the
rectangle specified should maintain an aspect ratio of 4x3, lest the
input video be distorted.

# Required arguments:

`input_video_filename`: char/string. Currently tested with 8-bit
grayscale .avi files and 24-bit color .mp4 files.

# Optional arguments (enter as name-value pairs after required
arguments):

`input_dist`: positive scalar in meters. Distance from lens to object(s)
in video.

`input_focal_length`: positive scalar in millimeters. Focal length of
video camera.

`input_sensor_size`: 1x2 positive scalar array in millimeters. Height
and width of video camera sensor.

`input_pixel_size`: positive scalar in microns. Height/width of square
pixels of video camera sensor.

`cropRect`: 1x4 nonnegative integer array in pixels. Top-left and
bottom-right corner pixels of rectangular area of input video to crop.
Aspect ratio should match `output_dim` or video will be distorted during
processing. [x1 y1 x2 y2].

`output_retina_size`: 1x2 positive scalar array in microns. Height and
width of target retina. Default [880 660] in Wei lab rig.

`output_viewing_dist`: positive scalar in centimeters. Optimal viewing
distance of animal to scene. Default 10.

`output_dim`: 1x2 positive integer array in pixels. Default [800 600] in
Wei lab OLED.

`output_resolution`: positive integer in pixels. Default 50. Size of
background flicker checkerboard squares. `output_dim` must be a multiple
of `output_resolution`.

`visualAngle2retinaDist`: positive scalar in microns per degree.
Conversion from visual angle to distance on retina. Default 32 for
mouse.

`input_timeframe`: 1x2 positive scalar array in either timeframes or
seconds. Defaults to entire duration of input video.

`output_dir`: char/string. Directory to write output videos. Default
input video file directory.

`kernel_factor`: scalar. Scaling factor applied to 3x3
convolution kernel: 
-1 -1 -1
-1  8 -1
-1 -1 -1
Default 2.

# Outputs:

Scaled video, convolved scaled video, and downsampled video. Written to
8-bit grayscale .avi files, as well as folders containing all video
frames. Default writes to same directory as input video file.
%}

function video2stim(input_filename, varargin)

%% Check inputs.
% Check input type and format. Assign default values if input parameter not provided.
p = inputParser;
v = @validateattributes;

addRequired( p,'input_filename',                    @(x) v(x,{'char','string'},{'nonempty'},mfilename,'input_filename',1));          % path to video file to process.
addParameter(p,'input_dist',                    [], @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'input_dist'));             % distance from camera lens to objects in scene (meters).
addParameter(p,'input_focal_length',            [], @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'input_focal_length'));     % camera focal length (millimeters).
addParameter(p,'input_sensor_size',             [], @(x) v(x,{'numeric'},{'numel',2,'positive'},mfilename,'input_sensor_size'));     % camera sensor size (millimeters).
addParameter(p,'input_pixel_size',              [], @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'input_pixel_size'));       % height/width of square pixels in camera sensor (microns).
addParameter(p,'cropRect',                      [], @(x) v(x,{'numeric'},{'numel',4,'nonnegative'},mfilename,'cropRect'));           % top-left and bottom-right pixel indices of rectangular area to crop. [x1 y1 x2 y2].
addParameter(p,'output_retina_size',     [880 660], @(x) v(x,{'numeric'},{'numel',2,'positive'},mfilename,'output_retina_size'));    % size of target retina (microns).
addParameter(p,'output_viewing_dist',           10, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'output_viewing_dist'));    % optimal viewing distance (centimeters).
addParameter(p,'output_dim',             [800 600], @(x) v(x,{'numeric'},{'numel',2,'positive'},mfilename,'output_dim'));            % height and width of output videos (pixels).
addParameter(p,'output_resolution',             50, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'output_resolution'));      % size of background flicker checkerboard squares (pixels).
addParameter(p,'visualAngle2retinaDist',        32, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'visualAngle2retinaDist')); % distance on target retina per visual angle (microns/degree).
addParameter(p,'input_timeframe',               [], @(x) v(x,{'numeric'},{'numel',2,'nonnegative'},mfilename,'input_timeframe'));    % temporal segment of input video to process (seconds or timeframes).
addParameter(p,'output_dir',                    [], @(x) v(x,{'char','string'},{'nonempty'},mfilename,'output_dir'));                % directory to write output videos.
addParameter(p,'kernel_factor',                  2, @(x) v(x,{'numeric'},{'scalar'},mfilename,'kernel_factor'));                     % multiplier to convolution kernel. Modulates overall intensity.

parse(p, input_filename, varargin{:});

% Write parsed inputs to variables.
input_dist = p.Results.input_dist;
input_focal_length = p.Results.input_focal_length;
input_sensor_size = p.Results.input_sensor_size;
input_pixel_size = p.Results.input_pixel_size;
cropRect = p.Results.cropRect;
output_retina_size = p.Results.output_retina_size;
output_viewing_dist = p.Results.output_viewing_dist;
output_dim = p.Results.output_dim;
output_resolution = p.Results.output_resolution;
visualAngle2retinaDist = p.Results.visualAngle2retinaDist;
input_timeframe = p.Results.input_timeframe;
output_dir = p.Results.output_dir;
kernel_factor = p.Results.kernel_factor;
clearvars p v varargin

% Check whether requested background flicker checkerboard square size can tile output video without gaps.
if any(floor(output_dim ./ output_resolution) ~= (output_dim ./ output_resolution))
    error('Pixel dimensions of output display/monitor must be a multiple of output resolution.')
end

% Check shape of array inputs.
if size(input_sensor_size,1) == 2
    input_sensor_size = input_sensor_size';
end
if size(output_retina_size,1) == 2
    output_retina_size = output_retina_size';
end
if size(output_dim,1) == 2
    output_dim = output_dim';
end
if size(cropRect,1) == 4
    cropRect = cropRect';
elseif size(cropRect,1) == 2
    cropRect(1,3:4) = cropRect(2,1:2);
    cropRect(2,:) = [];
end

% Convert to SI.
input_focal_length = input_focal_length * 10^-3;
input_sensor_size = input_sensor_size * 10^-3;
input_pixel_size = input_pixel_size * 10^-6;
output_retina_size = output_retina_size * 10^-6;
output_viewing_dist = output_viewing_dist * 10^-2;
visualAngle2retinaDist = visualAngle2retinaDist * 10^-6;

%% Prepare input video.
input_videoObj = VideoReader(input_filename);
clearvars input_filename

% Optionally, segment video temporally.
if isempty(input_timeframe)
    input_timeframe = [1 input_videoObj.NumFrames];
elseif ~isempty(input_timeframe)
    str = input('Is input timeframe in units of frame-numbers (y) or seconds (n)? y/n: ','s');
    while str ~= 'y' && str ~= 'n'
        str = input('Is input timeframe in units of frame-numbers (y) or seconds (n)? Please enter y/n: ','s');
    end
    if str == 'n'
        input_timeframe = input_timeframe .* input_fps + 1;
    end
    clearvars str
    if input_timeframe(2) > input_videoObj.NumFrames
        warning('Input video is not long enough to accomodate requested timeframe. Reading until end of input video.')
        input_timeframe(2) = input_videoObj.NumFrames;
    end
end

%% Determine pixel range to crop (optionally).
isCrop = 1;
if isempty(cropRect) && isempty(input_dist) && isempty(input_focal_length) && isempty(input_sensor_size) && isempty(input_pixel_size)
    warning('Processing entire frames without cropping.')
    isCrop = 0;
elseif ~isempty(cropRect)
    if cropRect(3) > input_videoObj.Width || cropRect(4) > input_videoObj.Height || cropRect(1) > input_videoObj.Width || cropRect(2) > input_videoObj.Height
        error('Area to crop specified by `cropRect` exceeds input video size %ix%i.', input_videoObj.Width, input_videoObj.Height)
    else
        % Check whether `cropRect` follows [topleft_x topleft_y botright_x botright_y].
        if cropRect(1) > cropRect(3) && cropRect(2) > cropRect(4)
            tmp = cropRect(1:2);
            cropRect(1:2) = cropRect(3:4);
            cropRect(3:4) = tmp;
            clearvars temp
        elseif cropRect(1) < cropRect(3) && cropRect(2) > cropRect(4)
            cropRect = [cropRect(1) cropRect(4) cropRect(3) cropRect(2)];
        elseif cropRect(1) > cropRect(3) && cropRect(2) < cropRect(4)
            cropRect = [cropRect(3) cropRect(2) cropRect(1) cropRect(4)];
        end
    end
else
    if isempty(input_dist) || isempty(input_focal_length) || (isempty(input_sensor_size) && isempty(input_pixel_size))
        error('Both `input_dist` and `input_focal_length`, as well as either `input_sensor_size` or `input_pixel_size`, is required to calculate area to crop using camera specifications.')
    elseif ~isempty(input_sensor_size) && ~isempty(input_pixel_size)
        if any([input_videoObj.Width input_videoObj.Height] * input_pixel_size ~= input_sensor_size)
            warning('Both `input_sensor_size` and `input_pixel_size` are given, but input video dimensions multiplied by `input_pixel_size` do not equal `input_sensor_size`. Using `input_pixel_size` to calculate area to crop.')
        end
    elseif isempty(input_pixel_size)
        input_pixel_size = input_sensor_size ./ [input_videoObj.Width input_videoObj.Height];
    end
    if isempty(output_retina_size) || isempty(visualAngle2retinaDist) || isempty(output_viewing_dist)
        error('`output_retina_size`, `visualAngle2retinaDist`, and `output_viewing_dist` are required to calculate area to crop based on camera specifications, target retina, and optimal viewing distance.')
    end
    clearvars input_sensor_size    
    
    % Find desired size of imagery depicted in output video.
    output_landscape_size = output_retina_size / visualAngle2retinaDist * pi/180 * output_viewing_dist; % x,y

    % Approximate size of imagery depicted in input video using pinhole model.
    input_landscape_per_pixel = input_pixel_size / input_focal_length * input_dist; % x,y

    % Find pixel dimensions to crop input video.
    crop_dim = round(output_landscape_size ./ input_landscape_per_pixel); % x,y.
    clearvars output_landscape_size input_landscape_per_pixel

    % Check if desired crop dimensions are larger than input video dimensions.
    if any(crop_dim > [input_videoObj.Width input_videoObj.Height])
        error('Landscape depicted in input video is not sufficiently large to accomodate the desired output landscape size determined by `output_retina_size`, `output_viewing_dist`, and `visualAngle2retinaDist`.');
    % If input video dimensions are larger than the desired cropped dimensions, query user for area to crop.
    elseif any(crop_dim < [input_videoObj.Width input_videoObj.Height])
        fprintf('Input video dimensions are larger than the desired cropped dimensions. \nInput video is %ix%i pixels, whereas the desired crop area is %ix%i pixels. \nChoose area to crop by dragging rectangle on pop-up window. Double-click rectangle when finished. \nNote that if the area of the rectangle is changed, area to crop may not correspond to requested input viewing distance or retina size. \nIf the aspect ratio of the rectangle is changed, video will be distorted during processing./n', input_videoObj.Width, input_videoObj.Height, crop_dim(1), crop_dim(2));
        testframe = squeeze(read(input_videoObj,1));
        f = figure; imshow(testframe);                                                 % display first frame of input video.
        r = drawrectangle('position',[0 0 crop_dim(1) crop_dim(2)],'Deletable',false); % overlay adjustable rectangle with dimensions of desired crop.
        wait(r);                                                                       % wait for user to double-click rectangle after dragging rectangle to select area to crop.
        close(f);
        cropRect = r.Position;
        clearvars r f testframe
        % Check that `cropRect` does not exceed video dimensions.
        cropRect(1:2) = round(cropRect(1:2));
        cropRect(1) = max(1,cropRect(1)); cropRect(2) = max(1,cropRect(2));
        if cropRect(1) + cropRect(3) > input_videoObj.Width
            cropRect(1) = cropRect(1) - 1;
        end
        if cropRect(2) + cropRect(4) > input_videoObj.Height
            cropRect(2) = cropRect(2) - 1;
        end
        cropRect(3:4) = cropRect(1:2) + cropRect(3:4);
    else
        isCrop = 0;
    end
    clearvars crop_dim
end
clearvars output_retina_size visualAngle2retinaDist input_focal_length input_dist input_pixel_size output_viewing_dist

%% Make directories to store processed video frames.
[~, filename] = fileparts(input_videoObj.Name);
mkdir(input_videoObj.Path, [filename '_crop_scale']);
mkdir(input_videoObj.Path, [filename '_crop_scale_filter']);
mkdir(input_videoObj.Path, [filename '_bkgdFlicker']);

%% Process video to produce frames.
isColor = strcmp(input_videoObj.VideoFormat, 'RGB24');
kernel = kernel_factor * ...
         [-1 -1 -1;
          -1  8 -1;
          -1 -1 -1]; % naive kernel for edge detection.
      
% Process frame-by-frame.
for i = input_timeframe(1):input_timeframe(2)
    frame = squeeze(read(input_videoObj,i)); % read frame from input video.
    if isCrop
        frame = frame(max(1,cropRect(1)-1):min(input_videoObj.Width,cropRect(3)+1), max(1,(cropRect(2)-1)):min(input_videoObj.Height,cropRect(4)+1),:); % crop frame leaving 1-pixel cushion to facilitate 3x3 kernel convolution.
    end
    if isColor
        frame = rgb2gray(frame); % convert to grayscale if necessary.
    end
    frame_filter = uint8(conv2(frame, kernel, 'same')); % apply high-pass filter.
    % If necessary, remove 1-pixel cushion.
    if isCrop
        if cropRect(1) > 1
            frame(:,1) = [];
            frame_filter(:,1) = [];
        end
        if cropRect(2) > 1
            frame(1,:) = [];
            frame_filter(1,:) = [];
        end
        if cropRect(3) < input_videoObj.Width
            frame(:,end) = [];
            frame_filter(:,end) = [];
        end
        if cropRect(4) < input_videoObj.Height
            frame(end,:) = [];
            frame_filter(end,:) = [];
        end
    end
    
    % Scale to match output pixel dimensions.
    frame = imresize(frame, [output_dim(2) output_dim(1)]);
    frame_filter = imresize(frame_filter, [output_dim(2) output_dim(1)]);
    
    % Down-sample to produce background checkerboard flicker.
    frame_downsample = imresize(frame, [output_dim(2)/output_resolution output_dim(1)/output_resolution]);
    bkgdFlicker = repelem(frame_downsample, output_resolution, output_resolution);
    
    % Write processed frames to disk.
    imgname = sprintf('frame%i.tif',i);
    imwrite(frame, fullfile(input_videoObj.Path, [filename '_crop_scale'], imgname), 'tif'); % write cropped and scaled frame.
    imwrite(frame_filter, fullfile(input_videoObj.Path, [filename '_crop_scale_filter'], imgname), 'tif'); % write cropped, scaled, and filtered frame.
    imwrite(bkgdFlicker, fullfile(input_videoObj.Path, [filename '_bkgdFlicker'], imgname), 'tif'); % write background flicker frame.
end
clearvars iscolor kernel frame frame_filter bkgdFlicker frame_downsample

%% Combine written frames into AVI video files.
if isempty(output_dir)
    output_dir = input_videoObj.Path;
end

outfile = VideoWriter(fullfile(output_dir, [filename '_crop_scale.avi']), 'Grayscale AVI');
outfile.FrameRate = input_videoObj.FrameRate;
open(outfile);

for i = input_timeframe(1):input_timeframe(2)
    frame = imread(fullfile(input_videoObj.Path, [filename '_crop_scale'], sprintf('frame%i.tif',i)));
    writeVideo(outfile, frame);
end
close(outfile);

outfile = VideoWriter(fullfile(output_dir, [filename '_crop_scale_filter.avi']), 'Grayscale AVI');
outfile.FrameRate = input_videoObj.FrameRate;
open(outfile);

for i = input_timeframe(1):input_timeframe(2)
    frame = imread(fullfile(input_videoObj.Path, [filename '_crop_scale_filter'], sprintf('frame%i.tif',i)));
    writeVideo(outfile, frame);
end
close(outfile);

outfile = VideoWriter(fullfile(output_dir, [filename '_bkgdFlicker.avi']), 'Grayscale AVI');
outfile.FrameRate = input_videoObj.FrameRate;
open(outfile);

for i = input_timeframe(1):input_timeframe(2)
    frame = imread(fullfile(input_videoObj.Path, [filename '_bkgdFlicker'], sprintf('frame%i.tif',i)));
    writeVideo(outfile, frame);
end
close(outfile);

end

