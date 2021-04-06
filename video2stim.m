%{
Convert input video to Wei lab rig OLED display stimulus. Match size of
depicted scene to size of target retina and optimal viewing distance.
Convolve with 3x3 kernel to detect edges. Scale video size to match OLED
size: 800x600 pixels. Downscale to match bipolar cell RF diamater to
produce background flicker.

# Required arguments:

`input_video_filename`: char. Currently tested with 8-bit grayscale .avi
files and 24-bit color .mp4 files.

`input_dist`: positive scalar in meters. Distance from lens to object(s)
in video.

`input_focal_length`: positive scalar in millimeters. Focal length of
video camera.

`input_sensor_size`: 1x2 positive scalar array in millimeters. Height
and width of video camera sensor.

# Optional arguments (enter as name-value pairs after required
arguments):

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

`output_dir`: char. Directory to write output videos. Default input
video file directory.

`kernel_factor`: optional scalar. Scaling factor applied to 3x3
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

function video2stim(input_filename, input_dist, input_focal_length, input_sensor_size, varargin)

%% Check inputs.
% Check input type and format. Assign default values if input parameter not provided.
p = inputParser;
v = @validateattributes;

addRequired( p,'input_filename',                    @(x) v(x,{'char','string'},{'nonempty'},mfilename,'input_filename',1));
addRequired( p,'input_dist',                        @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'input_dist',2));           % distance from camera lens to objects in scene (meters).
addRequired( p,'input_focal_length',                @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'input_focal_length',3));   % camera focal length (millimeters).
addRequired( p,'input_sensor_size',                 @(x) v(x,{'numeric'},{'numel',2,'positive'},mfilename,'input_sensor_size',4));   % camera sensor size (millimeters).
addParameter(p,'output_retina_size',     [880 660], @(x) v(x,{'numeric'},{'numel',2,'positive'},mfilename,'output_retina_size'));    % size of target retina (microns).
addParameter(p,'output_viewing_dist',           10, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'output_viewing_dist'));    % optimal viewing distance (centimeters).
addParameter(p,'output_dim',             [800 600], @(x) v(x,{'numeric'},{'numel',2,'positive'},mfilename,'output_dim'));            % height and width of output videos (pixels).
addParameter(p,'output_resolution',             50, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'output_resolution'));      % size of background flicker checkerboard squares (pixels).
addParameter(p,'visualAngle2retinaDist',        32, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'visualAngle2retinaDist')); % distance on target retina per visual angle (microns/degree).
addParameter(p,'input_timeframe',               [], @(x) v(x,{'numeric'},{'numel',2,'nonnegative'},mfilename,'input_timeframe'));    % temporal segment of input video to process (seconds or timeframes).
addParameter(p,'output_dir',                    [], @(x) v(x,{'char'},{'nonempty'},mfilename,'output_dir'));                         % directory to write output videos.
addParameter(p,'kernel_factor',                  2, @(x) v(x,{'numeric'},{'scalar'},mfilename,'kernel_factor'));                     % multiplier to convolution kernel. Modulates overall intensity.

parse(p, input_filename, input_dist, input_focal_length, input_sensor_size, varargin{:});

% Write parsed inputs to variables.
output_retina_size = p.Results.output_retina_size;
output_viewing_distance = p.Results.output_viewing_dist;
output_dim = p.Results.output_dim;
output_resolution = p.Results.output_resolution;
visualAngle2retinaDist = p.Results.visualAngle2retinaDist;
input_timeframe = p.Results.input_timeframe;
output_dir = p.Results.output_dir;
kernel_factor = p.Results.kernel_factor;
clearvars p v varargin

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

% Convert to SI.
input_focal_length = input_focal_length * 10^-3;
input_sensor_size = input_sensor_size * 10^-3;
output_retina_size = output_retina_size * 10^-6;
output_viewing_distance = output_viewing_distance * 10^-2;
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

%% Determine pixel range to crop.
% Find desired size of imagery depicted in output video.
output_landscape_size = output_retina_size / visualAngle2retinaDist * pi/180 * output_viewing_distance; % x,y

% Approximate size of imagery depicted in input video using pinhole model.
input_landscape_size = input_sensor_size / input_focal_length * input_dist; % x,y
input_landscape_per_pixel = input_landscape_size ./ [input_videoObj.Width input_videoObj.Height]; % x,y

% Find pixel dimensions to crop input video.
crop_dim = round(output_landscape_size ./ input_landscape_per_pixel); % x,y.
clearvars output_landscape_size input_landscape_size input_landscape_per_pixel output_retina_size angle2retina_dist input_sensor_size input_focal_length input_dist

% Check if desired crop dimensions are larger than input video dimensions.
if any(crop_dim > [input_videoObj.Width input_videoObj.Height])
    error(['Input video is not sufficiently large to accomodate the desired output viewing distance(s): ', sprintf('%.2f ', output_viewing_distance*10^2), 'cm.']);
% If input video dimensions are larger than the desired cropped dimensions, query user for area to crop.
elseif any(crop_dim < [input_videoObj.Width input_videoObj.Height])
    prompt = sprintf('Input video dimensions are larger than the desired cropped dimensions. \nInput video is %ix%i pixels, whereas the desired crop is %ix%i pixels.\nSelect top-left pixel index to crop by entering a two-membered numeric array (max [%i,%i]),\nor enter ''center'', ''top'', ''bottom'', ''left'', ''right'', ''topleft'', ''topright'', ''bottomleft'', or ''bottomright''. \nIf left empty, default crop center of input video:\n', input_videoObj.Width, input_videoObj.Height, crop_dim(1), crop_dim(2), input_videoObj.Width-crop_dim(1), input_videoObj.Height-crop_dim(2));
    topleft_coord = input(prompt);
    if isempty(topleft_coord)
        topleft_coord = round(([input_videoObj.Width input_videoObj.Height]-crop_dim)/2);
    else
        while ~strcmp(topleft_coord,'center') && ~strcmp(topleft_coord,'top') && ~strcmp(topleft_coord,'bottom') && ~strcmp(topleft_coord,'left') && ~strcmp(topleft_coord,'right') && ~strcmp(topleft_coord,'topright') && ~strcmp(topleft_coord,'topleft') && ~strcmp(topleft_coord,'bottomright') && ~strcmp(topleft_coord,'bottomleft') && ~(isnumeric(topleft_coord) && numel(topleft_coord) == 2 && all(topleft_coord + crop_dim <= [input_videoObj.Height, input_videoObj.Width]))
            if isnumeric(topleft_coord) && any(topleft_coord + crop_dim > [input_videoObj.Width input_videoObj.Height])
                topleft_coord = input(sprintf('Input top-left pixel index is too far right/down to accomodate desired cropped dimensions. Re-enter top-left pixel index (max [%i,%i]),\nor press enter to crop center of input video:\n', input_dim(1)-crop_dim(i,1), input_dim(2)-crop_dim(i,2)));
            else
                disp('Input not recognized./n')
                topleft_coord = input(prompt);
            end
        end
        if isnumeric(topleft_coord)
            elseif strcmpi(topleft_coord,'center')
                topleft_coord = round(([input_videoObj.Width input_videoObj.Height]-crop_dim)/2);
            elseif strcmpi(topleft_coord,'top')
                topleft_coord = [round((input_videoObj.Width-crop_dim(1))/2) 1];
            elseif strcmpi(topleft_coord,'bottom')
                topleft_coord = [round((input_videoObj.Width-crop_dim(1))/2) input_videoObj.Height-crop_dim(2)];
            elseif strcmpi(topleft_coord,'left')
                topleft_coord = [1 round((input_videoObj.Height-crop_dim(2))/2)];
            elseif strcmpi(topleft_coord,'right')
                topleft_coord = [input_videoObj.Width-crop_dim(1) round((input_videoObj.Height-crop_dim(2))/2)];
            elseif strcmpi(topleft_coord,'topleft')
                topleft_coord = [1 1];
            elseif strcmpi(topleft_coord,'topright')
                topleft_coord = [input_videoObj.Width-crop_dim(1) 1];
            elseif strcmpi(topleft_coord,'bottomleft')
                topleft_coord = [1 input_videoObj.Height-crop_dim(2)];
            elseif strcmpi(topleft_coord,'bottomright')
                topleft_coord = [input_videoObj.Width-crop_dim(1) input_videoObj.Height-crop_dim(2)];
        end
    end
else
    topleft_coord = [1 1];
end
clearvars prompt

%% Make directories to store processed video frames.
[~, filename] = fileparts(input_videoObj.Name);
mkdir(input_videoObj.Path, [filename '_crop_scale']);
mkdir(input_videoObj.Path, [filename '_crop_scale_filter']);
mkdir(input_videoObj.Path, [filename '_bkgdFlicker']);

%% Process video to produce frames.
iscolor = strcmp(input_videoObj.VideoFormat, 'RGB24');
kernel = kernel_factor * ...
         [-1 -1 -1;
          -1  8 -1;
          -1 -1 -1]; % naive kernel for edge detection.
      
% Process frame-by-frame.
for i = input_timeframe(1):input_timeframe(2)
    frame = squeeze(read(input_videoObj,i)); % read frame from input video.
    frame = frame(topleft_coord(2)-1:(topleft_coord(2)+crop_dim(2)), topleft_coord(1)-1:(topleft_coord(1)+crop_dim(1)),:); % crop frame leaving 1-pixel cushion to facilitate 3x3 kernel convolution.
    if iscolor
        frame = rgb2gray(frame); % convert to grayscale if necessary.
    end
    frame_filter = uint8(conv2(frame, kernel, 'same')); % apply high-pass filter.
    frame_filter = frame_filter(2:end-1,2:end-1); % remove 1-pixel cushion.
    frame = frame(2:end-1,2:end-1); % remove 1-pixel cushion.
    
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

