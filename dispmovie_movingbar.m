function dispmovie(varargin)

%% Check inputs.
if nargin == 0
    error('At least one input required: video filename or variable');
end

p = inputParser;
v = @validateattributes;

addParameter(p, 'video_filename', [],  @(x) v(x,{'char'},{'nonempty'},mfilename,'video_filename'));
addParameter(p, 'videomat', [], @(x) v(x,{'numeric'},{'nonempty','nonnan','nonnegative'},mfilename,'videomat'));
addParameter(p, 'video_fps', [], @(x) v(x,{'numeric'},{'scalar','nonnegative'},mfilename,'video_fps'));
addParameter(p, 'screenNumber', [], @(x) v(x,{'numeric'},{'scalar','integer','nonnegative'},mfilename,'screenNumber'));
addParameter(p, 'movingbar', 0, @(x) v(x,{'numeric'},{'scalar','binary'},mfilename,'movingbar'));
addParameter(p, 'nReps', 1, @(x) v(x,{'numeric'},{'scalar','positive','integer'},mfilename,'nReps'));
addParameter(p, 'nDirections', 4, @(x) v(x,{'numeric'},{'scalar','positive','integer'},mfilename,'nDirections'));
addParameter(p, 'bar_width', 100, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'bar_width'));
addParameter(p, 'bar_length', 800, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'bar_length'));
addParameter(p, 'bar_speed', 300, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'bar_speed'));
addParameter(p, 'bar_color', 1, @(x) v(x,{'numeric'},{'2d','>=',0,'<=',1},mfilename,'bar_color'));

parse(p, varargin{:});

video_filename = p.Results.video_filename;
videomat = p.Results.videomat;
video_fps = p.Results.video_fps;
screenNumber = p.Results.screenNumber;
movingbar = p.Results.movingbar;
nReps = p.Results.nReps;
nDirections = p.Results.nDirections;
bar_width = p.Results.bar_width;
bar_length = p.Results.bar_length;
bar_speed = p.Results.bar_speed;
bar_color = p.Results.bar_color;

clearvars varargin p v

%% Import movie.
if ~bitxor(isempty(video_filename), isempty(videomat))
    error('Precisely one video filename or video variable is required.')
end

if isempty(videomat)
    videoobj = VideoReader(video_filename);
    videomat = squeeze(read(videoobj));
    video_fps = videoobj.Framerate;
    clearvars video_filename
elseif isempty(video_fps)
    warning('Input video FPS not given. Assuming 60 FPS.')
    video_fps = 60;
end

videodim = size(videomat);
video_ifi = 1 / video_fps;

%% Optionally, create moving bar matrix.
if movingbar
    directions = repmat(0:(360/nDirections):(360-(360/nDirections)),nReps,1);
    for iRep = 1:nReps
        directions(iRep,:) = directions(iRep,randperm(nDirections));
    end
    disp(directions);
    
    barmat = ones(bar_width,bar_length+2,2);
    barmat(:,:,1) = barmat(:,:,1) * bar_color;
    barmat(:,2:end-1,2) = 0;
end

%% Initiate Psychtoolbox.
PsychDefaultSetup(2);
KbReleaseWait;

if isempty(screenNumber)
    screenNumber = max(Screen('Screens'));
end

black = BlackIndex(screenNumber);

try
    % Configure Psychtoolbox.
    Screen('Preference', 'SkipSyncTests', 0);
    oldVerbosityLevel = Screen('Preference', 'Verbosity', 0);
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 1);
    
    %% Open and configure window.
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
    [screenXpx, screenYpx] = Screen('WindowSize', window);
    screen_ifi = Screen('GetFlipInterval',window);
    topPriorityLevel = MaxPriority(window);
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    if screen_ifi > video_ifi
        warning('Input video FPS exceeds screen FPS by %.2f%%. May result in inaccurate frame timing.', (video_fps - 1/screen_ifi) / video_fps * 100)
    end
    
    if videodim(1) ~= screenYpx || videodim(2) ~= screenXpx
        warning('Input movie (%i,%i) and output screen (%i,%i) do not have the same dimensions.', videodim(2), videodim(1), screenXpx, screenYpx);
    end
    
    %% Convert movie frames to textures.
    frametex = NaN(videodim(end),1);
    for i = 1:videodim(end)
        frametex(i) = Screen('MakeTexture', window, videomat(:,:,i));
    end
    clearvars videomat
    
    if movingbar
        bartex = Screen('MakeTexture', window, barmat);
    end
    
    %% Show movie.
    % Configure window to prepare to show movie.
    HideCursor(window);
    KbReleaseWait;
    Priority(topPriorityLevel);
    vbl = Screen('Flip', window);
    
    % Display movie.
    tic;
    for i = 1:videodim(end)
        Screen('DrawTexture', window, frametex(i));
        vbl = Screen('Flip', window, vbl + 0.75*video_ifi);
        Screen('Close', frametex(i));
    end
    toc

    %% Restore original Psychtoolbox settings and close screen.
    Screen('Preference', 'Verbosity', oldVerbosityLevel);
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Priority(0);
    sca;
catch
    sca;
    psychrethrow(psychlasterror);
end

end