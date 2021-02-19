% Flip intensity values around average intensity for a given video
% matrix.

function outmat = reverse_brightness(inmat)

mean_intensity = mean(inmat(:));
outmat = 2*mean_intensity - inmat;

end