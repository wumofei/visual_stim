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

parse(p, varargin{:});

video_filename = p.Results.video_filename;
videomat = p.Results.videomat;
video_fps = p.Results.video_fps;
screenNumber = p.Results.screenNumber;
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
    window = PsychImaging('OpenWindow', screenNumber, black);
    [screenXpx, screenYpx] = Screen('WindowSize', window);
    screen_ifi = Screen('GetFlipInterval',window);
    topPriorityLevel = MaxPriority(window);
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    if screen_ifi > video_ifi
        warning('Screen FPS is lower than input movie FPS by %.2f%%. May result in inaccurate frame timing.', (video_fps - 1/screen_ifi) / video_fps * 100)
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
        if KbCheck
            break
        end
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