# RyanColourCode
Matlab code for generating LMS colour stimuli

What you need for the S-cone stimuli:

Files you need (or that Matlab loads and will need to be able to find):
•	The processed calibration data, including the spectra from the screen. These processed files should contain spectra that have been resampled. They need to be resampled because the sampling points at which the spectra come out from the Ocean Optics spectrometer are at weird locations, and we need them to match the sampling points of the cone fundamentals (sampled on and at every1 nm) so that the matrices are the same dimensions and can be multiplied. 
The spectra will depend on which display you are using (ie. Viewpixx or Propixx), and whether you are presenting stereo/dichoptic stimuli through the goggles (need measurements from both eyes) or binocularly (single measurement without the goggles). The files are located in /Calibration and are named, for eg:
Viewpixx_Processed_cal_data_2_4_2016.mat.

•	The Stockman and Sharpe cone fundamentals. These can be downloaded from here:
http://www.cvrl.org/cones.htm. I have downloaded the .csv file for the fundamentals sampled at 1 nm for both 2 and 10 degrees. I have then processed these in Matlab and re-saved them as .mat files: StockmanSharpe_2deg_cone_fundamentals_1nm.mat or StockmanSharpe_10deg_cone_fundamentals_1nm.mat and can be found in the folder /Colour. 

Code.
LMS2RGB_Vpixx.m (found in /Calibration)
This provides an RGB triplet for a specified LMS cone excitation, and needs both the RGB spectra and the cone fundamentals described above. It is called by any code that draws stimuli specified in LMS activations (e.g. those given below). This works for either the Viewpixx or Propixx. It is based on the Stanford vistadisp suite code cone2RGB.m and findMaxConeScale.m. It also requires the vistadisp-master suite to be in the Matlab path. Note that it is currently written for stereoscopic/dichoptic stimuli so it averages the L and R eye spectra (because they are both so similar). It would need to be adapted if you want to use binocular stimuli and spectra obtained without the goggles.

TryToMakeConeStimulus2.m (found in /Colour)
This gives a simple example of generating a Gaussian blob specified in LMS cone contrast, either negative or positive polarity. Plots images showing the appearance of the stimuli.

MID_Dots_test_colour_2.m (found in /Motion In Depth Colour)
Draws a simple oscillating random dot stereogram in colour, specified in LMS coordinates.

MID_Dots_find_isoluminance_2.m (found in /Motion In Depth Colour)
An important piece of code using heterochromatic flicker photometry to find the isoluminant point of 2 opposite sides of a dimension in LMS space.
This uses circular dots flickering between the 2 opposite poles and is updated with each button press. There are 6 trials: 3 for each of the 2 eyes (and the results are averaged per eye). Again, this will need to be adjusted if you’re using binocular stimuli. 


