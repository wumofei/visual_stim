function dispmovie(movie_filename, varargin)

%% Check inputs.
p = inputParser;
v = @validateattributes;

addRequired(p,'movie_filename',   @(x) v(x,{'char'},{'nonempty'},mfilename,'movie_filename',1));
addOptional(p,'screenNumber', [], @(x) v(x,{'numeric'},{'scalar','integer','nonnegative'}));

parse(p, movie_filename, varargin{:});

screenNumber = p.Results.screenNumber;
clearvars p v

%% Configure Psychtoolbox.
PsychDefaultSetup(2);
KbReleaseWait;

if isempty(screenNumber)
    screenNumber = max(Screen('Screens'));
end

black = BlackIndex(screenNumber);

%% Display movie.
try
    Screen('Preference', 'SkipSyncTests', 0);
    oldVerbosityLevel = Screen('Preference', 'Verbosity', 0);
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 1);
    window = PsychImaging('OpenWindow', screenNumber, black);
    movie = Screen('OpenMovie', window, movie_filename);
    topPriorityLevel = MaxPriority(window);
    HideCursor(window);
    Priority(topPriorityLevel);
    Screen('PlayMovie', movie, 1);

    while ~KbCheck
        tex = Screen('GetMovieImage', window, movie);
        if tex<=0
            break;
        end
    
        Screen('DrawTexture', window, tex);
        Screen('Flip', window);
        Screen('Close', tex);
    end

    Screen('PlayMovie', movie, 0);
    Screen('CloseMovie', movie);
    Screen('Preference', 'Verbosity', oldVerbosityLevel);
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Priority(0);
    sca;
catch
    sca;
    psychrethrow(psychlasterror);
end

end