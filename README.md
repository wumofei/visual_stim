# Natural Motion Stimulus

Process video file to create and display visual stimulus.

# Table of Contents
<!--ts-->
   * [Natural Motion Stimulus](#natural-motion-stimulus)
   * [Table of Contents](#table-of-contents)
   * [Dependencies](#dependencies)
   * [Usage](#usage)
   * [Bugs](#bugs)

<!-- Added by: mofei, at: Wed Jan 27 15:41:33 CST 2021 -->

<!--te-->

# Dependencies

* MATLAB R2019a or R2020b.
* Psychophysics Toolbox 3.0.17.
* Performance varies with graphics hardware.

# Usage

Call `video2stim.m` on video file (currently tested with 8-bit grayscale `.avi` files only) to produce a downsampled checkerboard flicker file and a convolved edge detection file. Call `dispmovie.m` on the edge detetion file to display directly. Call `dispmovie_movingbar.m` on the downsampled file to display a bright bar moving across background checkerboard flicker.

# Bugs

* Psychtoolbox may throw a synchronization failure error when calling `dispmovie.m` or `dispmovie_movingbar.m`. As far as I can tell, this is a hardware problem and not much can be done by changing code. Try rerunning a few times, or using better graphics hardware. You may skip Psychtoolbox synchronization tests by setting `Screen('Preference','SkipSyncTests',1)`, but this may result in inaccurate stimulus presentation.
