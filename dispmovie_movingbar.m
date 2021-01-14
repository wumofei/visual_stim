function dispmovie_movingbar(varargin)

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
addParameter(p, 'nReps', 1, @(x) v(x,{'numeric'},{'scalar','positive','integer'},mfilename,'nReps'));
addParameter(p, 'nDirs', 4, @(x) v(x,{'numeric'},{'scalar','positive','integer'},mfilename,'nDirections'));
addParameter(p, 'directions', [], @(x) v(x,{'numeric'},{'2d','nonempty'},mfilename,'directions'));
addParameter(p, 'shuffle_directions', 0, @(x) v(x,{'logical'},[],mfilename,'shuffle_directions'));
addParameter(p, 'bar_width', 100, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'bar_width'));
addParameter(p, 'bar_length', 800, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'bar_length'));
addParameter(p, 'bar_speed', 300, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'bar_speed'));
addParameter(p, 'bar_color', 1, @(x) v(x,{'numeric'},{'2d','>=',0,'<=',1},mfilename,'bar_color'));
addParameter(p, 'barDispBuffer', 10, @(x) v(x,{'numeric'},{'scalar','nonnegative'},mfilename,'barDispBuffer'));
addParameter(p, 'multisample', 4, @(x) v(x,{'numeric'},{'scalar','nonnegative','integer'},mfilename,'multisample'));

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
barDispBuffer = p.Results.barDispBuffer;
multisample = p.Results.multisample;

clearvars varargin p v

if ~bitxor(isempty(video_filename), isempty(videomat))
    error('Precisely one video filename or video variable is required.')
end

%% Import movie.
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
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, black, [], [], [], [], multisample);
    [screenXpx, screenYpx] = Screen('WindowSize', window);
    [xcenter, ycenter] = RectCenter(windowRect);
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
    
    %% Prepare moving bar texture.
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
   
    barmat = ones(bar_width+2,bar_length+2,2);
    barmat(:,:,1) = barmat(:,:,1) * bar_color;
    barmat([1 end],[1 end],2) = 0;
    bartex = Screen('MakeTexture', window, barmat);
    barDestRect = [0 0 bar_length+2 bar_width];
    clearvars barmat
    
    %% Calculate bar initial position(s).
    videoRect_cornerAngle = atand(videodim(1)/videodim(2));
    dist_center2bar_init = NaN(nTrials,1);
    
    % Try to remove for loop?
    for i = 1:nTrials
        direction = directions(i);
        if direction <= videoRect_cornerAngle || direction >= 360-videoRect_cornerAngle || (direction >= 180-videoRect_cornerAngle && direction <= 180+videoRect_cornerAngle)
            dist_center2bar_init(i) = videodim(2)/2/abs(cosd(direction)) + bar_length/2 + min([bar_width/2,videodim(1)/2,videodim(2)/2]).*abs(tand(direction));
        else
            dist_center2bar_init(i) = videodim(1)/2/abs(sind(direction)) + bar_length/2;
            if direction ~= 90 && direction ~= 270
                dist_center2bar_init(i) = dist_center2bar_init(i) + min([bar_width/2,videodim(1)/2,videodim(2)/2])./abs(tand(direction));
            end
        end
    end
    dist_center2bar_init = dist_center2bar_init + barDispBuffer;
    
    %% Show stimulus.
    % Configure window to prepare to show movie.
    HideCursor(window);
    KbReleaseWait;
    Priority(topPriorityLevel);
    
    for i = 1:nTrials
        direction = directions(i);
        dist_center2bar = dist_center2bar_init(i);
        Beeper('high');
        time0 = Screen('Flip', window);
        currentTime = time0;
        tic;
        while ~KbCheck && dist_center2bar >= -dist_center2bar_init(i)
            bar_pos = [xcenter+dist_center2bar*cosd(direction) ycenter-dist_center2bar*sind(direction)];
            frame2display = mod(round((currentTime-time0+screen_ifi)/video_ifi),videodim(end));
            barDestRect = CenterRectOnPoint(barDestRect, bar_pos(1), bar_pos(2));
            
            Screen('DrawTexture', window, frametex(frame2display));
            Screen('DrawTexture', window, bartex, [], barDestRect, -direction);
            currentTime = Screen('Flip', window, currentTime + 0.6*screen_ifi);
            
            dist_center2bar = dist_center2bar - bar_speed * screen_ifi;
        end
        toc
        Screen('FillRect', window, black);
        Screen('Flip', window);
        WaitSecs(1);
    end
    Screen('Close');

    %% Restore original Psychtoolbox settings and close screen.
    Screen('Preference', 'Verbosity', oldVerbosityLevel);
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Priority(0);
    Beeper('medium');
    sca;
catch
    sca;
    psychrethrow(psychlasterror);
end

end