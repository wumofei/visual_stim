# Natural Motion Stimulus

Process video file to create and display visual stimulus.

# Table of Contents
<!--ts-->
   * [Natural Motion Stimulus](#natural-motion-stimulus)
   * [Table of Contents](#table-of-contents)
   * [Dependencies](#dependencies)
   * [Usage](#usage)
      * [Example workflow](#example-workflow)
      * [video2stim.m](#video2stimm)
      * [Auxiliary functions](#auxiliary-functions)
   * [Bugs](#bugs)
   * [To-do](#to-do)

<!-- Added by: mofei, at: Mon May 17 14:33:19 CDT 2021 -->

<!--te-->

# Dependencies

* Tested on MATLAB R2019a and R2020b.
    * Image Processing Toolbox add-on is required for MATLAB R2020b.
* Psychophysics Toolbox 3.0.17.
* Performance varies with graphics hardware.

# Usage

Call `video2stim.m` on video file (currently tested with 8-bit grayscale `.avi` and 24-bit color `.mp4` files only) to produce a downsampled checkerboard flicker file and a convolved edge detection file. Call `dispmovie.m` on the edge detection file to display directly. Call `dispmovie_movingbar.m` on the downsampled file to display a bright bar moving across background checkerboard flicker.

## Example workflow

```
>> cd visual_stim
>> video2stim('myvideo.avi', 'input_dist', 0.5, 'input_focal_length', 16, 'input_sensor_size', [36 24]);
% corresponds to 0.5 meter distance between camera and object, 16mm lens, and 36x24mm sensor size for standard full-frame 35mm cameras.
% Select area to crop by adjusting rectangular indicator in pop-up window. Double-click rectangle to confirm.
...
>> dispmovie('myvideo_crop_scale_filter.avi');
...
>> dispmovie_movingbar('myvideo_bkgdFlicker.avi');
...
```

## `video2stim.m`

Two ways of selecting region of input video to crop are supported:

1. Specify `input_dist`, `input_focal_length`, and either`input_sensor_size` or `input_pixel_size`. Along with `visualAngle2retinaDist`, `output_retina_size` and `output_viewing_dist`, pixel dimensions of input video that match requested optimal viewing distance are calculated in lines 219-226. A pop-up window displaying the first frame of the input video, overlaid with a rectangular indicator with the calculated crop pixel dimensions, will appear. Move rectangular indicator to select desired region to crop. Double-click to finish. You may change the size of the rectangular indicator, but keep in mind that you might want to keep the aspect ratio constant, lest the input video be distorted.

2. Directly specify a rectangular region in pixels to crop using `cropRect`. Keep in mind that the Wei lab OLED is 800x600 pixels, so the rectangle specified should maintain an aspect ratio of 4x3, lest the input video be distorted.

## Auxiliary functions

* `readmovie.m`: load video file into memory. Useful for repeated operations on small videos. Will exceed memory limits for large videos.
* `writemovie.m`: write video matrix from memory to disk.
* `reverse_brightness.m`: flip brightness around mean intensity for a given video matrix.
* `dispmovie_long.m`: show video without preloading entire video into memory. May result in more frame timing delays. Use if input video exceeds memory availability.
* `dispmovie_movingbar_long.m`: show moving bar over background video without preloading entire video into memory.
* `video_linstretch.m`: adjust mean luminance and contrast of video. Unfinished! Algorithm used to determine pixel value limits needs improvement, specifically to accommodate case where user wishes to change mean luminance but leave contrast unchanged at value `1`. Current algorithm involves dividing by `1 - contrast` on lines 69 and 98, which is `0` when `contrast = 1`. [Michelson contrast](https://en.wikipedia.org/wiki/Contrast_(vision)#Michelson_contrast) is used where $I_{min}$ is usually `0` for many videos, and thereby `contrast = 1`. See similar inbuilt MATLAB functions `imadjust.m` and `stretchlim.m` that also include step of determining pixel limits for possible ways to fix.

# Bugs

* Psychtoolbox may throw a synchronization failure error when calling `dispmovie.m` or `dispmovie_movingbar.m`. As far as I can tell, this is a hardware problem and not much can be done by changing code. Try rerunning a few times, or using better graphics hardware. You may skip Psychtoolbox synchronization tests by setting `Screen('Preference','SkipSyncTests',1)`, but this may result in inaccurate stimulus presentation.

* `video_linstretch.m` does not behave as wanted when desired new contrast value equals 1. Need to fix algorithm for determining pixel value limits.

# To-do

* Fix `video_linstretch.m`.
* Research edge detection kernels.
* Add recording triggers.
* Add gamma table lookup?
