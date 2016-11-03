

% Written in an attempt to work out the procedure to draw Gaussian blob stimulus
% for something specified in LMS cone activations, using the cal data for the Viewpixx, obtained 2/4/16
% LMS triplets have a max of 1.0 (not sure if min is -1 or 0).
% This relies on the function LMS2RGB_Viewpixx, which is based on the vista-disp code, cone2RGB & findMaxConeScale
% It follows a failed/aborted attempt to draw the stimuli in DKL space.

%R Maloney, April 2016


%Stimulus details from MID_Cue_thresholds_2IFC_5.m / standard MID dot stimulus
PPD = 37; 
% specify dot size:
dot_sigma_in_degrees = 0.1;            % 0.064; %size of SD of the dot profile in degs/vis angle
dot_sigma = dot_sigma_in_degrees * PPD; % sigma in pixels
dotsize = round(dot_sigma * 10);        % make the dots some multiple of sigma
% NOTE: dotsize is simply half the length of the sides of a square patch of pixels that the dot profile is placed within.
% It is dot_sigma that really determines the size of the dots. Obviously, dotsize must be larger than dot_sigma.

% make the dot profile:
x = (1-dotsize)/2:(dotsize-1)/2;
%x = x/PPD; %rescale x into vis angle
[x,y] = meshgrid(x,x);
y = -y;
[a,r] = cart2pol(x,y);

% This gives us a white-peaked dot (+ve contrast polarity)
env = [];
env = exp(-r.^2/(2*dot_sigma.^2)); % gaussian window
env = env./max(max(abs(env)));     % normalize peak to +/- 1
env2 = -env;  

img = repmat(env,1,1,3); %to pre-allocate the image matrix

% Load the processed cal data for the Viewpixx:
load('\\storage.its.york.ac.uk\pshome\rm1380\My Documents\Calibration\Viewpixx_Processed_cal_data_2_4_2016.mat')

% OR, load the processed cal data for the PROpixx: 
load ('\\PSHome\Home\rm1380\My Documents\Calibration\PROpixx_Processed_cal_data_21_4_2016.mat')

% Load the Stockman/Sharpe 2 deg cone fundamentals:
load('\\storage.its.york.ac.uk\pshome\rm1380\My Documents\Colour\StockmanSharpe_2deg_cone_fundamentals_1nm.mat')

theta = pi/2 + deg2rad(1); %theta should be somewhere around pi/2
%Compute RGB values for every pixel in LMS cone activation space.

stimLMS.dir = [cos(theta)/sqrt(2), cos(theta)/sqrt(2), sin(theta)] %the desired cone activation
%stimLMS.dir = [1 1 1] 
% LMS2RGB_Vpixx will set the LMS.scale to the max possible if you request something outside the gamut, so just set it to 1.0 if you want max contrast.
DesiredScale = 1 %the scale factor (akin to desired cone contrast, I think); max for s-cones is about 0.748 on the viewpixx
tic
for ii = 1:dotsize
    for jj = 1:dotsize
        
        stimLMS.scale = env(ii,jj) * DesiredScale; %Work out cone activation for this pixel
        stimRGB = LMS2RGB_Vpixx (stimLMS, fundamentals, resampledSpectra);
        img(ii, jj, :) = stimRGB.dir * stimRGB.scale + 0.5; %scale & add the background to the obtained RGB triplet
    end
end
toc

figure
imshow(img);
%put the background in.
%Note that currently, we simply use the default background of 0.5, 0.5, 0.5
background = [0.5 0.5 0.5]; 
set(gcf, 'color', background);







%-------------------- older method --------------------
%Remember that the 3rd dimension of env is for the L, M & S inputs, in that order!
% We multiply L by the first row of LMS2RGB
% envRGB(:,:,1) = env(:,:,1)*LMS2RGB(1) + env(:,:,1)*LMS2RGB(4) + env(:,:,1)*LMS2RGB(7);
% % And multiply M by the 2nd row of LMS2RGB
% envRGB(:,:,2) = env(:,:,2)*LMS2RGB(2) + env(:,:,2)*LMS2RGB(5) + env(:,:,2)*LMS2RGB(8);
% % And finally, multiply S by the 3rd row of LMS2RGB
% envRGB(:,:,3) = env(:,:,3)*LMS2RGB(3) + env(:,:,3)*LMS2RGB(6) + env(:,:,3)*LMS2RGB(9);
% 
% %Now, add 0.5 to all matrices so they modulate around  mid-grey:
% envRGB = envRGB + 0.5;
% 
% %(Normalise the peak of envRGB to 1.0, because 'image' wont plot for values outside 0,1
% envRGB(:,:,3) = envRGB(:,:,3)./max(max(abs(envRGB(:,:,3))));     % normalize peak to +/- 1
% 
% figure
% image(envRGB)
% axis off
% axis square