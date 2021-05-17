% Flip intensity values around average intensity for a given 8-bit
% grayscale video.

function outmat = reverse_brightness(input_video, varargin)

%% Check inputs.
% Check input type and format. Assign default values if input parameter not provided.
p = inputParser;
v = @validateattributes;

addRequired( p, 'input_video',     @(x) v(x,{'char','string','numeric'},{'nonempty'},mfilename,'input_video'));    % filename or matrix.
addParameter(p, 'video_fps',   [], @(x) v(x,{'numeric'},{'scalar','nonnegative'},mfilename,'video_fps')); % default 60 fps.

parse(p, input_video, varargin{:});

video_fps = p.Results.video_fps;
clearvars varargin p v

% If video filename is given, read video file. Otherwise, use input video matrix.
if ischar(input_video)
    videoObj = VideoReader(input_video);
    clearvars input_video
    if ~strcmp(videoObj.VideoFormat, 'Grayscale')
        error('Input video must be grayscale.') % assume function is called on post-processed video from video2stim.m
    end
    % Make directory to store processed frames.
    [~, filename] = fileparts(input_videoObj.Name);
    mkdir(videoObj.Path, [filename '_reverseBrightness']);
    % Process frame-by-frame.
    for i = 1:videoObj.NumFrames
        frame = squeeze(read(videoObj,i)); % read frame from input video.
        mean_intensity = mean(frame(:));
        frame = 2*mean_intensity - frame;
        
        % Write processed frame to disk.
        imgname = sprintf('frame%i.tif',i);
        imwrite(frame, fullfile(videoObj.Path, [filename '_reverseBrightness'], imgname), 'tif'); % write brightness reversed frame.
    end
    % Combine written frames into output movie.
    outfile = VideoWriter(fullfile(videoObj.Path, [filename '_reverseBrightness']), 'Grayscale AVI');
    if isempty(video_fps)
        outfile.FrameRate = videoObj.FrameRate;
    else
        outfile.FrameRate = video_fps;
    end
    open(outfile);
    for i = 1:videoObj.NumFrames
        frame = imread(fullfile(videoObj.Path, [filename '_reverseBrightness'], sprintf('frame%i.tif',i)));
        writeVideo(outfile,frame);
    end
    close(outfile);
else
    if ndims(input_video) ~= 3
        error('Input video must be grayscale.')
    end
    mean_intensity = mean(input_video(:));
    outmat = 2*mean_intensity - input_video;
end

end