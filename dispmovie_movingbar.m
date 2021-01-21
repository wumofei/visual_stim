% This function displays a single bar of uniform intensity moving at
% uniform speed and various directions across an arbitrary
% black-and-white movie background.

function dispmovie_movingbar(input_video, varargin)

%% Check inputs.
% Check input type and format. Assign default values if input parameter not provided.
p = inputParser;
v = @validateattributes;

addRequired( p, 'input_video',             @(x) v(x,{'char','numeric'},{'nonempty'},mfilename,'input_video'));
addParameter(p, 'video_filename',      [], @(x) v(x,{'char'},{'nonempty'},mfilename,'video_filename'));
addParameter(p, 'videomat',            [], @(x) v(x,{'numeric'},{'nonempty','nonnan','nonnegative'},mfilename,'videomat'));
addParameter(p, 'video_fps',           [], @(x) v(x,{'numeric'},{'scalar','nonnegative'},mfilename,'video_fps'));
addParameter(p, 'screenNumber',        [], @(x) v(x,{'numeric'},{'scalar','integer','nonnegative'},mfilename,'screenNumber'));
addParameter(p, 'nReps',                1, @(x) v(x,{'numeric'},{'scalar','positive','integer'},mfilename,'nReps'));
addParameter(p, 'nDirs',                4, @(x) v(x,{'numeric'},{'scalar','positive','integer'},mfilename,'nDirections'));
addParameter(p, 'directions',          [], @(x) v(x,{'numeric'},{'2d','nonempty'},mfilename,'directions'));
addParameter(p, 'shuffle_directions',   0, @(x) v(x,{'logical'},[],mfilename,'shuffle_directions'));
addParameter(p, 'bar_width',          100, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'bar_width'));
addParameter(p, 'bar_length',         800, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'bar_length'));
addParameter(p, 'bar_speed',          300, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'bar_speed'));
addParameter(p, 'bar_color',            1, @(x) v(x,{'numeric'},{'2d','>=',0,'<=',1},mfilename,'bar_color'));
addParameter(p, 'multisample',          8, @(x) v(x,{'numeric'},{'scalar','nonnegative','integer'},mfilename,'multisample'));

parse(p, varargin{:});

video_filename = p.Results.video_filename;
videomat = p.Results.videomat;
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

% Check whether unique input video is given.
if ~bitxor(isempty(video_filename), isempty(videomat))
    error('Precisely one video filename or video variable is required.')
end

%% Import movie.
% If video filename is given, read video file. Otherwise, use input video matrix.
if isempty(videomat)
    videoobj = VideoReader(video_filename);
    videomat = squeeze(read(videoobj));
    video_fps = videoobj.Framerate;
    clearvars video_filename
elseif isempty(video_fps)
    warning('Input video FPS not given. Assuming 60 FPS.');
    video_fps = 60;
end

videodim = size(videomat); % videodim(1) = height, videodim(2) = width, videodim(3) = timeframes. Assuming grayscale video.
video_ifi = 1 / video_fps; % interframe interval in seconds.

%% Initiate Psychtoolbox.
PsychDefaultSetup(2); % AssertOpenGL, unify key names, unify color ranges to 0 thru 1 (instead of 0 thru 255)
KbReleaseWait;

if isempty(screenNumber)
    screenNumber = max(Screen('Screens')); % find screen number corresponding to external screen. If connected to multiple external screens, find screen with largest port number.
end

black = BlackIndex(screenNumber); % usually equal to 0.

% Use try/catch to avoid getting stuck on a black screen.
try
    % Configure Psychtoolbox.
    Screen('Preference', 'SkipSyncTests', 0);                          % perform standard screen sync tests.
    oldVerbosityLevel = Screen('Preference', 'Verbosity', 0);          % suppress non-critical errors. Save old setting to variable.
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 1); % suppress loading screen. Save old setting to variable.
    
    %% Open and configure window.
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, black, [], [], [], [], multisample); % open black window, with multisample option to (hopefully) improve diagonal movement.
    [screenXpx, screenYpx] = Screen('WindowSize', window); % get screen size.
    [xcenter, ycenter] = RectCenter(windowRect); % get screen center.
    screen_ifi = Screen('GetFlipInterval',window); % interframe interval in seconds.
    topPriorityLevel = MaxPriority(window); % prepare CPU priority during stimulus presentation.
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA'); % standard configuration for blending.
    clearvars windowRect
    
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
    % Standard procedure to display a matrix of luminance values using Psychtoolbox is to:
    % (1) convert matrix to texture using Screen('MakeTexture')
    % (2) draw texture using Screen('DrawTexture')
    % (3) display on screen using Screen('Flip')
    % We proceed with step (1) here.
    frametex = NaN(videodim(end),1);
    for i = 1:videodim(end)
        frametex(i) = Screen('MakeTexture', window, videomat(:,:,i));
    end
    clearvars videomat
    
    %% Create circular mask.
    maskmat = zeros(screenYpx, screenXpx, 2); % maskmat(:,:,1) contains grayscale intensity values from 0 thru 1. maskmat(:,:,2) contains transparency values from 0 (transparent) thru 1 (opaque).
    maskradius = min([videodim(1); videodim(2); screenXpx; screenYpx])/2;
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
        dist_center2bar = dist_center2bar_init;           % write initial bar position to variable that updates during each while loop.
        
        Beeper('high');                 % make sound at beginning of each stim presentation trial.
        time0 = Screen('Flip', window); % record stim presentation start time.
        currentTime = time0;            % write initial time to variable that updates during each while loop.
        
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
            % screen_ifi. If moving bar presentation is longer than length of input
            % video, loop video.
            frame2display = mod(round((currentTime-time0+screen_ifi)/video_ifi),videodim(end));
            
            % Place moving bar destination on calculated moving bar center position.
            barDestRect = CenterRectOnPoint(barDestRect, bar_pos(1), bar_pos(2));
            
            % Draw textures, sequentially from background to foreground.
            Screen('DrawTexture', window, frametex(frame2display));
            Screen('DrawTexture', window, bartex, [], barDestRect, -direction);
            Screen('DrawTexture', window, masktex);
            
            % Flip to screen with option to ensure precise timing. Indicates all
            % drawing commands should be finished and ready to flip to screen by 0.7
            % of an interframe interval.
            currentTime = Screen('Flip', window, currentTime + 0.7*screen_ifi);
            
            % Update bar position.
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