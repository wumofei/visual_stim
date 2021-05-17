%{
Display movie using Psychtoolbox without preloading entire video into
memory, at the cost of possible frame timing error. Recommend trying
`dispmovie.m` first. If memory problems arise, then use
`dispmovie_long.m`.

# Required arguments: 

`input_video`: char. Filename of input video. Input video constrained to
grayscale by Wei lab OLED.

# Optional arguments (enter as name-value pairs after required
arguments):

`video_fps`: scalar in frames per second. Default 60. To present sped up
or slowed down video, adjust `video_fps` accordingly. Keep in mind that
Wei lab OLED usually operates at 60Hz refresh rate. Ideally `video_fps`
would be equal to a multiple of 60, or vice versa, to ensure relatively
accurate frame presentation timing.

`screenNumber`: integer. Number corresponding to screen on which to
display stimulus. Default finds largest screennumber among connected
screens.
%}

function dispmovie_long(input_video, varargin)

%% Check inputs.
% Check input type and format. Assign default values if input parameter not provided.
p = inputParser;
v = @validateattributes;

addRequired( p, 'input_video',      @(x) v(x,{'char','string'},{'nonempty'},mfilename,'input_video'));                           % filename.
addParameter(p, 'video_fps',    [], @(x) v(x,{'numeric'},{'scalar','nonnegative'},mfilename,'video_fps'));              % default 60 fps.
addParameter(p, 'screenNumber', [], @(x) v(x,{'numeric'},{'scalar','integer','nonnegative'},mfilename,'screenNumber')); % default find external screen with greatest port number.

parse(p, input_video, varargin{:});

video_fps = p.Results.video_fps;
screenNumber = p.Results.screenNumber;
clearvars varargin p v

%% Import movie.
videoObj = VideoReader(input_video);
if ~strcmp(videoObj.VideoFormat, 'Grayscale')
    error('Input video must be grayscale.')
end
if isempty(video_fps)
    video_fps = videoObj.Framerate;
end
clearvars input_video

video_ifi = 1 / video_fps; % interframe interval in seconds.

%% Initiate Psychtoolbox.
PsychDefaultSetup(2); % AssertOpenGL, unify key names, set color ranges to 0 thru 1 (instead of 0 thru 255).
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
    window = PsychImaging('OpenWindow', screenNumber, black);                  % open blank screen.
    [screenXpx, screenYpx] = Screen('WindowSize', window);                     % get screen size.
    screen_ifi = Screen('GetFlipInterval',window);                             % interframe interval in seconds.
    topPriorityLevel = MaxPriority(window);                                    % prepare CPU priority during stimulus presentation.
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA'); % standard configuration for blending.
    
    % By default, input video is centered on the OLED screen, whose
    % dimensions are 800x600 pixels. Please scale/crop input video
    % appropriately before calling this function.
    if videoObj.Height ~= screenYpx || videoObj.Width ~= screenXpx
        warning('Input movie (%ix%i) and output screen (%ix%i) do not have the same dimensions.', videoObj.Width, videoObj.Height, screenXpx, screenYpx);
    end
    
    %% Show movie.
    % Configure window to prepare to show movie.
    HideCursor(window);
    KbReleaseWait;
    Priority(topPriorityLevel);
    
    % Show first frame.
    frame = squeeze(read(videoObj,1));
    frametex = Screen('MakeTexture', window, frame);
    Screen('DrawTexture', window, frametex);
    time0 = Screen('Flip', window); % mark stim presentation start time.
    currentTime = time0; % write initial time to variable that updates during each while loop.
    
    % Display movie.
    while ~KbCheck && currentTime - time0 + screen_ifi <= videoObj.NumFrames * video_ifi
        % Find index of video frame to display based on estimate of next screen
        % update. Elapsed time since start of stimulus is given by
        % currentTime-time0. Best estimate of time until next screen flip is
        % screen_ifi. Frame to display is frame corresponding to midpoint of
        % next screen interval, i.e. elapsed time + 1.5 * screen interval.
        frame2display = floor((currentTime-time0+1.5*screen_ifi)/video_ifi)+1;
        
        frame = squeeze(read(videoObj, frame2display));
        frametex = Screen('MakeTexture', window, frame);
        Screen('DrawTexture', window, frametex);
        
        % Flip to screen with option to ensure precise timing. Indicates all
        % drawing commands should be finished and ready to flip to screen by 0.7
        % of an interframe interval.
        currentTime = Screen('Flip', window, currentTime + 0.7*video_ifi);
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