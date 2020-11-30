function dispmovie(movie_filename, varargin)

%% Check inputs.
p = inputParser;
v = @validateattributes;

addRequired(p,'movie_filename',   @(x) v(x,{'char'},{'nonempty'},mfilename,'movie_filename',1));
addOptional(p,'screenNumber', [], @(x) v(x,{'numeric'},{'scalar','integer','nonnegative'}));

parse(p, movie_filename, varargin{:});

screenNumber = p.Results.screenNumber;
clearvars p v

%% Prepare display.
PsychDefaultSetup(2);
KbReleaseWait;

if nargin < 2 || isempty(screenNumber)
    screenNumber = max(Screen('Screens'));
end

try
    window = PsychImaging('OpenWindow', screenNumber, 0);
    movie = Screen('OpenMovie', window, movie_filename);
    topPriorityLevel = MaxPriority(window);
    HideCursor(window);

    %% Display movie.
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
    Priority(0);
    sca;
catch
    sca;
    psychrethrow(psychlasterror);
end

end