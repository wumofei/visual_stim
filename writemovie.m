% Write matrix in memory to video file. Default format grayscale avi.

function writemovie(videomat, fps, filename, format)

if nargin < 4
    outfile = VideoWriter(filename, 'Grayscale AVI');
else
    outfile = VideoWriter(filename, format);
end

outfile.FrameRate = fps;
open(outfile);
writeVideo(outfile, videomat);
close(outfile);

end