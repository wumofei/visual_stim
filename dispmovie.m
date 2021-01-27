%{
Display movie using Psychtoolbox.

# Required arguments: 

`input_video`: char or matrix. If char, interpret as video filename.
Otherwise, interpret as 8-bit grayscale video matrix. Currently tested
with 8-bit grayscale .avi files only.

# Optional arguments (enter as name-value pairs after required
arguments):

`video_fps`: scalar in frames per second. Default 60.

`screenNumber`: integer. Number corresponding to screen on which to
display stimulus. Default finds largest screennumber among connected
screens.
%}

function dispmovie(input_video, varargin)

%% Check inputs.
% Check input type and format. Assign default values if input parameter not provided.
p = inputParser;
v = @validateattributes;

addRequired( p, 'input_video',      @(x) v(x,{'char','numeric'},{'nonempty'},mfilename,'input_video'));                 % filename or matrix.
addParameter(p, 'video_fps',    [], @(x) v(x,{'numeric'},{'scalar','nonnegative'},mfilename,'video_fps'));              % default 60 fps.
addParameter(p, 'screenNumber', [], @(x) v(x,{'numeric'},{'scalar','integer','nonnegative'},mfilename,'screenNumber')); % % default find external screen with greatest port number.

parse(p, input_video, varargin{:});

video_fps = p.Results.video_fps;
screenNumber = p.Results.screenNumber;
clearvars varargin p v

%% Import movie.
% If video filename is given, read video file. Otherwise, use input video matrix.
if ischar(input_video)
    videoobj = VideoReader(input_video);
    videomat = squeeze(read(videoobj));
    video_fps = videoobj.Framerate;
else
    videomat = input_video;
    if isempty(video_fps)
        warning('Input video FPS not given. Assuming 60 FPS.');
        video_fps = 60;
    end
end
clearvars input_video

videodim = size(videomat); % videodim(1) = height, videodim(2) = width, videodim(3) = timeframes. Assuming grayscale video.
video_ifi = 1 / video_fps; % interframe interval in seconds.

%% Initiate Psychtoolbox.
PsychDefaultSetup(2); % AssertOpenGL, unify key names, unify color ranges to 0 thru 1 (instead of 0 thru 255).
KbReleaseWait;

if isempty(screenNumber)
    screenNumber = max(Screen('Screens')); % find screen number corresponding to external screen. If connected to multiple external screens, find screen with largest port number.
end

black = BlackIndex(screenNumber); % usually equal to 0.

% Use try/catch to avoid getting stuck on a black screen upon error.
try
    % Configure Psychtoolbox.
    Screen('Preference', 'SkipSyncTests', 0);                          % perform standard screen sync tests.
    oldVerbosityLevel = Screen('Preference', 'Verbosity', 0);          % suppress non-critical errors. Save old setting to variable.
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 1); % suppress loading screen. Save old setting to variable.
    
    %% Open and configure window.
    window = PsychImaging('OpenWindow', screenNumber, black); % open blank screen.
    [screenXpx, screenYpx] = Screen('WindowSize', window);                     % get screen size.
    screen_ifi = Screen('GetFlipInterval',window);                             % interframe interval in seconds.
    topPriorityLevel = MaxPriority(window);                                    % prepare CPU priority during stimulus presentation.
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA'); % standard configuration for blending.
    
    % Wei lab OLED usually operates at a refresh rate of 60 Hz. If input
    % video FPS is greater than 60, some frames will be skipped.
    if screen_ifi > video_ifi
        warning('Input video FPS exceeds screen refresh rate by %.2f%%. May result in inaccurate frame timing.', (video_fps - 1/screen_ifi) / video_fps * 100)
    end
    
    % By default, input video is centered on the OLED screen, whose
    % dimensions are 800x600 pixels. Please scale/crop input video
    % appropriately before calling this function.
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
    vbl = Screen('Flip', window); % mark stim presentation start time.
    
    % Display movie.
    for i = 1:videodim(end)
        Screen('DrawTexture', window, frametex(i));
        
        % Flip to screen with option to ensure precise timing. Indicates all
        % drawing commands should be finished and ready to flip to screen by 0.7
        % of an interframe interval.
        vbl = Screen('Flip', window, vbl + 0.7*video_ifi);
        
        if KbCheck
            break
        end
    end
    Screen('Close');

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