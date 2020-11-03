%{
Parameters:

input_video_filename: char.

input_distance: positive scalar in units of meters.

input_focal_length: positive scalar in units of millimeters.

input_sensor_size: 1x2 positive scalar array in units of millimeters.

output_retinal_size: 1x2 positive scalar array in units of microns. Default [880 660].

output_viewing_distance: 1xn positive scalar array in units of centimeters. Default 10.

output_dim: 1x2 positive integer array in units of pixels. Default [800 600].

output_resolution: positive integer in units of pixels. Default 50.

visual_angle2retinal_length: positive scalar in units of microns per degree. Default 32.

input_timeframe: optional 1x2 positive scalar array in units of either
time-frames or seconds. If left empty, defaults to entire duration of
input video.

output_directory: optional char. If left empty, defaults to current directory.

kernel_factor: optional scalar. Scaling factor applied to 3x3 convolution kernel
-1 -1 -1
-1  8 -1
-1 -1 -1

Outputs:

Scaled video, convolved scaled video, and downsampled video.Written to .avi files for each viewing distance given.
%}

function video2stim(input_video_filename, input_distance, input_focal_length, input_sensor_size, output_retinal_size, output_viewing_distance, output_dim, output_resolution, visual_angle2retinal_length, input_timeframe, output_directory, kernel_factor)

%% Check inputs.
% Check number of inputs and input type.
if nargin < 4
    error('At least 4 input arguments required.')
elseif ~ischar(input_video_filename)
    error('Input video path must be entered as a char.')
elseif ~isscalar(input_distance)
    error('Input video distance must be entered as a scalar (in units of meters).')
elseif ~isscalar(input_focal_length)
    error('Input video focal length must be entered as a scalar (in units of millimeters).')
elseif (numel(input_sensor_size) ~= 2) || (~isnumeric(input_sensor_size))
    error('Input video sensor size (in both x and y dimensions) must be entered as a two-membered numeric array (in units of millimeters).')
end

if nargin < 5 || isempty(output_retinal_size)
    output_retinal_size = [880 660];
elseif numel(output_retinal_size) ~= 2 || ~isnumeric(output_retinal_size)
    error('Output retinal size (in both x and y dimensions) must be entered as a two-membered numeric array (in units of microns on target retina).')
end

if nargin < 6 || isempty(output_viewing_distance)
    output_viewing_distance = 10;
elseif ~isnumeric(output_viewing_distance)
    error('Output optimal viewing distance must be numeric (in units of centimeters).')
end

if nargin < 7 || isempty(output_dim)
    output_dim = [800 600];
elseif numel(output_dim) ~= 2 || ~isnumeric(output_dim)
    error('Output dimensions (in both x and y dimensions) must be entered as a two-membered numeric array (in units of pixel).')
end

if nargin < 8 || isempty(output_resolution)
    output_resolution = 50;
elseif ~isscalar(output_resolution)
    error('Output resolution must be entered as a scalar (in units of pixels).')
end

if nargin < 9 || isempty(visual_angle2retinal_length)
    visual_angle2retinal_length = 32;
elseif ~isscalar(visual_angle2retinal_length)
    error('Conversion between one degree visual angle and microns on target retina must be entered as a scalar.')
end

if nargin >= 10 && ~isempty(input_timeframe)
    if numel(input_timeframe) ~= 2 || ~isnumeric(input_timeframe)
        error('Input timeframe must be entered as a two-membered numeric array (in units of seconds or frame-numbers).')
    end
end

if nargin < 11 || isempty(output_directory)
    output_directory = pwd;
elseif ~ischar(output_directory)
    error('Output directory must be entered as a char.')
end

if nargin < 12 || isempty(kernel_factor)
    kernel_factor = 2;
elseif ~isscalar(kernel_factor) || kernel_factor <= 0
    error('High-pass filter kernel factor must be a positive scalar.')
end

% Check shape of array parameters.
if size(input_sensor_size,1) == 2
    input_sensor_size = input_sensor_size';
end
if size(output_retinal_size,1) == 2
    output_retinal_size = output_retinal_size';
end
if size(output_viewing_distance,2) ~= 1
    output_viewing_distance = output_viewing_distance';
end
if size(output_dim,1) == 2
    output_dim = output_dim';
end

% Convert to SI.
input_focal_length = input_focal_length * 10^-3;
input_sensor_size = input_sensor_size * 10^-3;
output_retinal_size = output_retinal_size * 10^-6;
output_viewing_distance = output_viewing_distance * 10^-2;
visual_angle2retinal_length = visual_angle2retinal_length * 10^-6;

%% Read input video.
input_videoobj = VideoReader(input_video_filename);
[~,filename,~] = fileparts(input_video_filename);
input_fps = input_videoobj.Framerate;
clear input_video_filename

% Optionally, segment video temporally.
if nargin < 10 || isempty(input_timeframe)
    input_timeframe = [1 input_videoobj.Duration*input_fps]; % select full-length input video
else
    str = input('Is input timeframe in units of frame-numbers (y) or seconds (n)? y/n: ','s');
    while str ~= 'y' && str ~= 'n'
        str = input('Is input timeframe in units of frame-numbers (y) or seconds (n)? Please enter y/n: ','s');
    end
    if str == 'n'
        input_timeframe = input_timeframe .* input_fps;
    end
end
clear str

% Extract frames from video object.
% input_videoobj.CurrentTime = (input_timeframe - 1) / input_videoobj.Framerate;
input_video = squeeze(read(input_videoobj, input_timeframe));
input_dim = [size(input_video,2) size(input_video,1)]; % x,y
nFrames = size(input_video,3);
clear input_timeframe input_videoobj

%% Apply high-pass filter.
kernel = kernel_factor * ...
         [-1 -1 -1;
          -1  8 -1;
          -1 -1 -1];
input_video_filtered = zeros(input_dim(2),input_dim(1),nFrames, 'uint8');

for i = 1:nFrames
    input_video_filtered(:,:,i) = imfilter(input_video(:,:,i), kernel, 'conv');
end
clear kernel

%% Crop input video to obtain desired output video size.
% Find desired size of imagery depicted in output video.
output_landscape_size = output_retinal_size / visual_angle2retinal_length * pi/180 .* output_viewing_distance; % x,y

% Approximate size of imagery depicted in input video using pinhole model.
input_landscape_size = input_sensor_size / input_focal_length * input_distance; % x,y
input_landscape_per_pixel = input_landscape_size ./ input_dim; % x,y

% Find pixel dimensions to crop input video.
crop_dim = output_landscape_size ./ input_landscape_per_pixel; % x,y. Non-integer.
crop_dim = round(crop_dim);
clear output_landscape_size input_landscape_size input_landscape_per_pixel output_retinal_size visual_angle2retinal_length input_sensor_size input_focal_length input_distance

%{
% Make desired crop dimensions integers that are scaled appropriately to the desired output pixel dimensions.
reduced_output_dim = output_dim / gcd(output_dim(1),output_dim(2));
rounded_crop_dim = reduced_output_dim .* round(crop_dim ./ reduced_output_dim);
for i = 1:size(crop_dim,1)
   if abs(rounded_crop_dim(i,1)-crop_dim(i,1)) < abs(rounded_crop_dim(i,2)-crop_dim(i,2))
        crop_dim(i,1) = rounded_crop_dim(i,1);
        crop_dim(i,2) = rounded_crop_dim(i,1) / reduced_output_dim(1) * reduced_output_dim(2);
    else
        crop_dim(i,2) = rounded_crop_dim(i,2);
        crop_dim(i,1) = rounded_crop_dim(i,2) / reduced_output_dim(2) * reduced_output_dim(1);
    end
end
%}

% Check if desired crop dimensions are larger than input video dimensions.
if any(crop_dim > input_dim)
    idInsufficient = sum(crop_dim > input_dim,2) > 0;
    disp(['Input video is not sufficiently large to accomodate the following desired output viewing distance(s): ', sprintf('%.2f ', output_viewing_distance(idInsufficient)*10^2), 'cm.']);
    crop_dim(idInsufficient) = [];
end
if isempty(crop_dim)
    error('Input video is not sufficiently large to accomodate any desired viewing distance.')
end

nOutputs = size(crop_dim,1);
clear idInsufficient

% If input video dimensions are larger than the desired cropped dimensions, query user for area to crop.
topleft_coord = ones(nOutputs,2);
if any(crop_dim < input_dim)
    for i = 1:nOutputs
        if sum(crop_dim(i,:) == input_dim,2) ~= 2
            prompt = sprintf('Input video dimensions are larger than the desired cropped dimensions corresponding to desired viewing distance %.2f centimeters.\nInput video is %ix%i pixels, whereas the desired crop is %ix%i pixels.\nSelect top-left pixel index to crop by entering a two-membered numeric array (max [%i,%i]),\nor enter ''center'', ''top'', ''bottom'', ''left'', ''right'', ''topleft'', ''topright'', ''bottomleft'', or ''bottomright''. \nIf left empty, default crop center of input video:\n', output_viewing_distance(i)*10^2, input_dim(1), input_dim(2), crop_dim(i,1), crop_dim(i,2), input_dim(1)-crop_dim(i,1), input_dim(2)-crop_dim(i,2));
            topleft_coord_i = input(prompt);
            if isempty(topleft_coord_i)
                topleft_coord_i = round((input_dim-crop_dim(i,:))/2);
            else
                while ~strcmp(topleft_coord_i,'center') && ~strcmp(topleft_coord_i,'top') && ~strcmp(topleft_coord_i,'bottom') && ~strcmp(topleft_coord_i,'left') && ~strcmp(topleft_coord_i,'right') && ~strcmp(topleft_coord_i,'topright') && ~strcmp(topleft_coord_i,'topleft') && ~strcmp(topleft_coord_i,'bottomright') && ~strcmp(topleft_coord_i,'bottomleft') && ~(isnumeric(topleft_coord_i) && numel(topleft_coord_i) == 2 && all(topleft_coord_i + crop_dim(i,:) <= input_dim))
                    if isnumeric(topleft_coord_i) && any(topleft_coord_i + crop_dim(i,:) > input_dim)
                        topleft_coord_i = input(sprintf('Input top-left pixel index is too far right/down to accomodate desired cropped dimensions. Re-enter top-left pixel index (max [%i,%i]),\nor press enter to crop center of input video:\n', input_dim(1)-crop_dim(i,1), input_dim(2)-crop_dim(i,2)));
                    else
                        disp('Input not recognized./n')
                        topleft_coord_i = input(prompt);
                    end
                end
                if isnumeric(topleft_coord_i)
                elseif strcmp(topleft_coord_i,'center')
                    topleft_coord_i = round((input_dim-crop_dim(i,:))/2);
                elseif strcmp(topleft_coord_i,'top')
                    topleft_coord_i = [round((input_dim(1)-crop_dim(i,1))/2) 1];
                elseif strcmp(topleft_coord_i,'bottom')
                    topleft_coord_i = [round((input_dim(1)-crop_dim(i,1))/2) input_dim(2)-crop_dim(2)];
                elseif strcmp(topleft_coord_i,'left')
                    topleft_coord_i = [1 round((input_dim(2)-crop_dim(i,2))/2)];
                elseif strcmp(topleft_coord_i,'right')
                    topleft_coord_i = [input_dim(1)-crop_dim(1) round((input_dim(2)-crop_dim(i,2))/2)];
                elseif strcmp(topleft_coord_i,'topleft')
                    topleft_coord_i = [1 1];
                elseif strcmp(topleft_coord_i,'topright')
                    topleft_coord_i = [input_dim(1)-crop_dim(1) 1];
                elseif strcmp(topleft_coord_i,'bottomleft')
                    topleft_coord_i = [1 input_dim(2)-crop_dim(2)];
                elseif strcmp(topleft_coord_i,'bottomright')
                    topleft_coord_i = [input_dim(1)-crop_dim(1) input_dim(2)-crop_dim(2)];
                end
            end
        end
        topleft_coord(i,:) = topleft_coord_i;
    end
end
clear topleft_coord_i prompt input_dim

% Crop input video.
cropped_video = zeros(max(crop_dim(:,2)), max(crop_dim(:,1)), nFrames, nOutputs, 'uint8');
cropped_video_filtered = cropped_video;
% cropped_video_filtered_normalized = cropped_video;
for i = 1:nOutputs
    cropped_video(1:crop_dim(i,2),1:crop_dim(i,1),:,i) = input_video(topleft_coord(2):(topleft_coord(2)+crop_dim(i,2)-1), topleft_coord(1):(topleft_coord(1)+crop_dim(i,1)-1), :);
    cropped_video_filtered(1:crop_dim(i,2),1:crop_dim(i,1),:,i) = input_video_filtered(topleft_coord(2):(topleft_coord(2)+crop_dim(i,2)-1), topleft_coord(1):(topleft_coord(1)+crop_dim(i,1)-1), :);
%     for j = 1:nFrames
%         cropped_video_filtered_normalized(:,:,j,i) = histeq(cropped_video_filtered(:,:,j,i), imhist(cropped_video(:,:,j,i)));
%     end
end

%{
for i = 1:nOutputs
    output_file = VideoWriter([output_directory '/' filename '_out_cropped_viewingdistance' num2str(output_viewing_distance(i)*10^2,2) 'cm.avi'], 'Grayscale AVI');
    output_file.FrameRate = input_fps;
    open(output_file);
    writeVideo(output_file, cropped_video(:,:,:,i));
    close(output_file);
end
%}

clear topleft_coord input_video_filtered

%% Scale cropped video to obtain pixel dimensions of output display/monitor.
scaled_video = zeros(output_dim(2), output_dim(1), nFrames, nOutputs, 'uint8');
scaled_video_filtered = scaled_video;
% scaled_video_filtered_normalized = scaled_video;
for i = 1:nOutputs
    scaled_video(:,:,:,i) = imresize(cropped_video(1:crop_dim(i,2), 1:crop_dim(i,1),:,i), [output_dim(2) output_dim(1)]);
    scaled_video_filtered(:,:,:,i) = imresize(cropped_video_filtered(1:crop_dim(i,2), 1:crop_dim(i,1),:,i), [output_dim(2) output_dim(1)]);
%     scaled_video_filtered_normalized(:,:,:,i) = imresize(cropped_video_filtered_normalized(1:crop_dim(i,2), 1:crop_dim(i,1),:,i), [output_dim(2) output_dim(1)]);
end

% Write output videos.
for i = 1:nOutputs
    output_file = VideoWriter([output_directory '/' filename '_out_scaled_viewingdistance' num2str(output_viewing_distance(i)*10^2,2) 'cm.avi'], 'Grayscale AVI');
    output_file.FrameRate = input_fps;
    open(output_file);
    writeVideo(output_file, scaled_video(:,:,:,i));
    close(output_file);
    
    output_file = VideoWriter([output_directory '/' filename '_out_scaled_filtered_viewingdistance' num2str(output_viewing_distance(i)*10^2,2) 'cm.avi'], 'Grayscale AVI');
    output_file.FrameRate = input_fps;
    open(output_file);
    writeVideo(output_file, scaled_video_filtered(:,:,:,i));
    close(output_file);
    
%     output_file = VideoWriter([output_directory '/' filename '_out_scaled_filtered_normalized_viewingdistance' num2str(output_viewing_distance(i)*10^2,2) 'cm.avi'], 'Grayscale AVI');
%     output_file.FrameRate = input_fps;
%     open(output_file);
%     writeVideo(output_file, scaled_video_filtered_normalized(:,:,:,i));
%     close(output_file);
end
clear crop_dim cropped_video nFrames cropped_video_filtered cropped_video_filtered_normalized


%% Down-sample to obtain desired (reduced) output resolution.
downsample_video = imresize(scaled_video, [output_dim(2)/output_resolution output_dim(1)/output_resolution]);
% downsample_video_filtered = imresize(scaled_video_filtered, [output_dim(2)/output_resolution output_dim(1)/output_resolution]);
% downsample_video_filtered_normalized = imresize(scaled_video_filtered_normalized, [output_dim(2)/output_resolution output_dim(1)/output_resolution]);
output_video = repelem(downsample_video, output_resolution, output_resolution);
% output_video_filtered = repelem(downsample_video_filtered, output_resolution, output_resolution);
% output_video_filtered_normalized = repelem(downsample_video_filtered_normalized, output_resolution, output_resolution);
clear output_dim downsample_video output_resolution scaled_video scaled_video_filtered downsample_video_filtered scaled_video_filtered_normalized downsample_video_filtered_normalized

%% Write output video to .avi file.    
for i = 1:nOutputs
    output_file = VideoWriter([output_directory '/' filename '_out_final_viewingdistance' num2str(output_viewing_distance(i)*10^2,2) 'cm.avi'], 'Grayscale AVI');
    output_file.FrameRate = input_fps;
    open(output_file);
    writeVideo(output_file, output_video(:,:,:,i));
    close(output_file);
    
%     output_file = VideoWriter([output_directory '/' filename '_out_final_filtered_viewingdistance' num2str(output_viewing_distance(i)*10^2,2) 'cm.avi'], 'Grayscale AVI');
%     output_file.FrameRate = input_fps;
%     open(output_file);
%     writeVideo(output_file, output_video_filtered(:,:,:,i));
%     close(output_file);
    
%     output_file = VideoWriter([output_directory '/' filename '_out_final_filtered_normalized_viewingdistance' num2str(output_viewing_distance(i)*10^2,2) 'cm.avi'], 'Grayscale AVI');
%     output_file.FrameRate = input_fps;
%     open(output_file);
%     writeVideo(output_file, output_video_filtered_normalized(:,:,:,i));
%     close(output_file);
end

end

