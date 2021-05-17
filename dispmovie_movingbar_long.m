%{
Display a single bar of uniform intensity moving at uniform speed
along various directions across an arbitrary black-and-white movie
background with a circular mask.

Does not preload video into memory, at the cost of possible frame timing
error. Recommend trying `dispmovie_movingbar.m` first. If memory
problems arise, then use `dispmovie_movingbar_long.m`.

# Required arguments:

`input_video`: char. Filename of background video. Input video
constrained to grayscale by Wei lab OLED.

# Optional arguments (enter as name-value pairs after required
arguments):

`video_fps`: scalar in frames per second. Default 60. To present sped up
or slowed down video, adjust `video_fps` accordingly. Keep in mind that
Wei lab OLED usually operates at 60Hz refresh rate. Ideally `video_fps`
would be a multiple of 60, or vice versa, to ensure relatively accurate
frame presentation timing.

`screenNumber`: integer. Number corresponding to screen on which to
display stimulus. Default finds largest screennumber among connected
screens.

`nReps`: integer. Number of repetitions to display stimulus. Default 1.

`nDirs`: integer. During each repetition, number of directions to show
moving bar. Default 4.

`directions`: numeric 2D array. Specify directions to show moving bar.
If left empty, finds `nDirs` number of directions equally spaced from 0
to 360 degrees. Default empty.

`shuffle_directions`: binary scalar. Permute directions during each
repetition. Default 0.

`bar_width`: positive scalar in pixels. Default 100.

`bar_length`: positive scalar in pixels. Default 800.

`bar_speed`: positive scalar in pixels per second. Default 300.

`bar_color`: scalar between 0 and 1. Grayscale intensity of moving bar.

`multisample`: positive integer. Reduce jagged artifacts during diagonal
bar movement by calculating multiple samples for each pixels.
Performance depends on graphics hardware. Default 4.

%}
function dispmovie_movingbar_long(input_video, varargin)

%% Check inputs.
% Check input type and format. Assign default values if input parameter not provided.
p = inputParser;
v = @validateattributes;

addRequired( p, 'input_video',             @(x) v(x,{'char','string'},{'nonempty'},mfilename,'input_video'));                           % filename of background video.
addParameter(p, 'video_fps',           [], @(x) v(x,{'numeric'},{'scalar','nonnegative'},mfilename,'video_fps'));              % default 60 fps.
addParameter(p, 'screenNumber',        [], @(x) v(x,{'numeric'},{'scalar','integer','nonnegative'},mfilename,'screenNumber')); % default find external screen with greatest port number.
addParameter(p, 'nReps',                1, @(x) v(x,{'numeric'},{'scalar','positive','integer'},mfilename,'nReps'));           % repetitions of stimulus to show.
addParameter(p, 'nDirs',                4, @(x) v(x,{'numeric'},{'scalar','positive','integer'},mfilename,'nDirections'));     % within each repetition, directions to show. E.g., `nDirs = 4` means shoing bar moving at 0, 90, 180, and 270 degrees.
addParameter(p, 'directions',          [], @(x) v(x,{'numeric'},{'2d','nonempty'},mfilename,'directions'));                    % optionally specify directions of moving bar.
addParameter(p, 'shuffle_directions',   0, @(x) v(x,{'logical'},[],mfilename,'shuffle_directions'));                           % whether to randomly permute directions within each repetition.
addParameter(p, 'bar_width',          100, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'bar_width'));                 % pixels.
addParameter(p, 'bar_length',         800, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'bar_length'));                % pixels.
addParameter(p, 'bar_speed',          300, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'bar_speed'));                 % pixels per second.
addParameter(p, 'bar_color',            1, @(x) v(x,{'numeric'},{'2d','>=',0,'<=',1},mfilename,'bar_color'));                  % grayscale intensity of moving bar. Between 0 and 1.
addParameter(p, 'multisample',          4, @(x) v(x,{'numeric'},{'scalar','nonnegative','integer'},mfilename,'multisample'));  % samples to calculate per pixel to reduce jagged diagonal movement. Depends on graphics hardware.

parse(p, input_video, varargin{:});

% Write parsed inputs.
video_fps = p.Results.video_fps;
screenNumber = p.Results.screenNumber;
nReps = p.Results.nReps;
nDirs = p.Results.nDirs;
directions = p.Results.directions;
shuffle_directions = p.Results.shuffle_directions;
bar_width = p.Results.bar_width;
bar_length = p.Results.bar_length;
bar_speed = p.Results.bar_speed;
bar_color = p.Results.bar_color;
multisample = p.Results.multisample;

clearvars varargin p v

%% Import movie.
videoObj = VideoReader(input_video);
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
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, black, [], [], [], [], multisample); % open black window, with multisample option to (hopefully) improve diagonal movement.
    [screenXpx, screenYpx] = Screen('WindowSize', window);                     % get screen size.
    [xcenter, ycenter] = RectCenter(windowRect);                               % get screen center.
    screen_ifi = Screen('GetFlipInterval',window);                             % interframe interval in seconds.
    topPriorityLevel = MaxPriority(window);                                    % prepare CPU priority during stimulus presentation.
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA'); % standard configuration for blending.
    clearvars windowRect
    
    % By default, input video is centered on the OLED screen, whose
    % dimensions are 800x600 pixels. Please scale/crop input video
    % appropriately before calling this function.
    if videoObj.Height ~= screenYpx || videoObj.Width ~= screenXpx
        warning('Input movie (%ix%i) and output screen (%ix%i) do not have the same dimensions.', videoObj.Width, videoObj.Height, screenXpx, screenYpx);
    end
    
    %% Create circular mask.
    maskmat = zeros(screenYpx, screenXpx, 2); % maskmat(:,:,1) contains grayscale intensity values from 0 thru 1. maskmat(:,:,2) contains transparency values from 0 (transparent) thru 1 (opaque).
    maskradius = min([videoObj.Height; videoObj.Width; screenXpx; screenYpx])/2;
    for x = 1:screenXpx
        for y = 1:screenYpx
            if (x-xcenter)^2+(y-ycenter)^2 > maskradius^2
                maskmat(y,x,2) = 1;
            end
        end
    end
    masktex = Screen('MakeTexture', window, maskmat);
    clearvars screenXpx screenYpx maskmat
    
    %% Prepare moving bar texture.
    barmat = ones(bar_width+2,bar_length+2,2);      % barmat(:,:,2) again contains transparency values.
    barmat(:,:,1) = barmat(:,:,1) * bar_color;
    barmat([1 end],[1 end],2) = 0;                  % set edge transparency to 0 (completely transparent) to facilitate smoother motion.
    bartex = Screen('MakeTexture', window, barmat);
    barDestRect = [0 0 bar_length+2 bar_width];     % initiate destination position to draw moving bar.
    clearvars barmat
    
    %% Prepare moving bar directions.
    if isempty(directions)
        directions = repmat(0:(360/nDirs):(360-(360/nDirs)),nReps,1);
    end
    nDirs = length(directions);
    if shuffle_directions
        for i = 1:nReps
            directions(i,:) = directions(i,randperm(nDirs));
        end
    end
    disp('Directions: ');
    disp(directions);
    directions = directions(:);
    nTrials = length(directions);
    
    %% Show stimulus.
    % Configure window to prepare to show movie.
    HideCursor(window);
    KbReleaseWait;
    Priority(topPriorityLevel);
    
    for i = 1:nTrials
        direction = directions(i);
        direction_cos = cosd(direction);
        direction_sin = sind(direction);
        
        % Start moving bar just outside circular mask.
        dist_center2bar_init = maskradius + bar_length/2; % distance from screen center to starting bar center.
        dist_center2bar = dist_center2bar_init; % write initial bar position to variable that updates during each while loop.
        
        Beeper('high'); % make sound at beginning of each stim presentation trial.
        
        % Show initial frame.
        bar_pos = [xcenter+dist_center2bar*direction_cos, ycenter-dist_center2bar*direction_sin]; % Find position of moving bar center.
        barDestRect = CenterRectOnPoint(barDestRect, bar_pos(1), bar_pos(2)); % Place moving bar destination on calculated moving bar center position.
        frame = squeeze(read(videoObj,1)); % Read first frame of input video.
        frametex = Screen('MakeTexture', window, frame); % Convert frame to texture.
        Screen('DrawTexture', window, frametex); % Show first frame of input video.
        Screen('DrawTexture', window, bartex, [], barDestRect, -direction); % Draw moving bar.
        Screen('DrawTexture', window, masktex); % Draw circular mask.
        time0 = Screen('Flip', window); % record stim presentation start time.
        currentTime = time0; % write initial time to variable that updates during each while loop.
        dist_center2bar = dist_center2bar - bar_speed * screen_ifi; % Update bar distance from center of screen.
        
        % Stim presentation loop.
        % Break when keyboard pressed. End presentation when bar completely exits circular mask.
        while ~KbCheck && dist_center2bar >= -dist_center2bar_init
            
            % Find position of moving bar center. Since indexing along y-axis is top
            % to bottom (e.g. top-left corner is [0,0], bottom-left corner is [0,
            % 600]), subtract from ycenter.
            bar_pos = [xcenter+dist_center2bar*direction_cos, ycenter-dist_center2bar*direction_sin];
            
            % Find index of video frame to display based on estimate of next screen
            % update. Elapsed time since start of stimulus is given by
            % currentTime-time0. Best estimate of time until next screen flip is
            % screen_ifi. Frame to display is frame corresponding to midpoint of
            % next screen interval, i.e. elapsed time + 1.5 * screen interval. If
            % moving bar presentation is longer than length of input video, loop
            % video using mod.
            frame2display = mod(floor((currentTime-time0+1.5*screen_ifi)/video_ifi)+1,videoObj.NumFrames);
            
            % Place moving bar destination on calculated moving bar center position.
            barDestRect = CenterRectOnPoint(barDestRect, bar_pos(1), bar_pos(2));
            
            % Draw textures, sequentially from background to foreground.
            frame = squeeze(read(videoObj, frame2display));
            frametex = Screen('MakeTexture', window, frame);
            Screen('DrawTexture', window, frametex);
            Screen('DrawTexture', window, bartex, [], barDestRect, -direction);
            Screen('DrawTexture', window, masktex);
            
            % Flip to screen with option to ensure precise timing. Indicates all
            % drawing commands should be finished and ready to flip to screen by 0.7
            % of an interframe interval. Step (3) of process mentioned above.
            currentTime = Screen('Flip', window, currentTime + 0.7*screen_ifi);
            
            % Update bar distance from center of screen.
            dist_center2bar = dist_center2bar - bar_speed * screen_ifi;
        end
        
        % Pause between different direction trials.
        Screen('FillRect', window, black);
        Screen('Flip', window);
        WaitSecs(1);
    end
    Screen('Close');

    %% Restore original Psychtoolbox settings and close screen.
    Screen('Preference', 'Verbosity', oldVerbosityLevel);
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Priority(0);
    Beeper('medium'); % make sound at end of stim presentation.
    sca;
catch
    sca;
    psychrethrow(psychlasterror);
end

end