# Natural Motion Stimulus

Process video file to create and display visual stimulus.

# Table of Contents
<!--ts-->
   * [Natural Motion Stimulus](#natural-motion-stimulus)
   * [Table of Contents](#table-of-contents)
   * [Dependencies](#dependencies)
   * [Usage](#usage)
      * [Example workflow](#example-workflow)
      * [Auxiliary functions](#auxiliary-functions)
   * [Bugs](#bugs)
   * [To-do](#to-do)

<!-- Added by: mofei, at: Fri Feb 26 11:47:30 CST 2021 -->

<!--te-->

# Dependencies

* Tested on MATLAB R2019a and R2020b. Image Processing Toolbox add-on is required for MATLAB R2020b.
* Psychophysics Toolbox 3.0.17.
* Performance varies with graphics hardware.

# Usage

Call `video2stim.m` on video file (currently tested with 8-bit grayscale `.avi` and 24-bit color `.mp4` files only) to produce a downsampled checkerboard flicker file and a convolved edge detection file. Call `dispmovie.m` on the edge detection file to display directly. Call `dispmovie_movingbar.m` on the downsampled file to display a bright bar moving across background checkerboard flicker.

## Example workflow

```
>> cd visual_stim
>> video2stim('myvideo.avi', 0.5, 16, [36 24]); % corresponds to 0.5 meter distance between camera and object, 16mm lens, and 36x24mm sensor size for standard full-frame 35mm cameras.
...
>> dispmovie('myvideo_crop_scale_filter.avi');
...
>> dispmovie_movingbar('myvideo_bkgdFlicker.avi');
...
```

## Auxiliary functions

* `readmovie.m`: load video file into memory. Useful for repeated operations on small videos. Will exceed memory limits for large video.
* `writemovie.m`: write video matrix from memory to disk.
* `reverse_brightness.m`: flip brightness around mean intensity for a given video matrix.

# Bugs

* Psychtoolbox may throw a synchronization failure error when calling `dispmovie.m` or `dispmovie_movingbar.m`. As far as I can tell, this is a hardware problem and not much can be done by changing code. Try rerunning a few times, or using better graphics hardware. You may skip Psychtoolbox synchronization tests by setting `Screen('Preference','SkipSyncTests',1)`, but this may result in inaccurate stimulus presentation.

* May run into memory problems on large video files.

# To-do

* Reduce memory usage.
* Integrate mean intensity and contrast adjustments in `video_linstretch.m` into display functions.
* Research edge detection kernels.
* Add recording triggers.
* Add gamma table lookup?
