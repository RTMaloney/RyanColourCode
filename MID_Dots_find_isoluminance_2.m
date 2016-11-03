
function MID_Dots_find_isoluminance_2 (SubjCode, runN)

% Works on either the PROpixx or the Viewpixx.

% Version _2: done using DrawDots to save time & memory.
% Ie not drawing the Gaussian blobs, since each pixel needs to be allocated RGB coordinates:
% this was taking about 1 sec per dot to run through LMS2RGB (ie both colours).

%%%%---------------------------------%%%%
%     Set some initial parameters:
%%%%---------------------------------%%%%

UseS_cone = true; %Flag for S cone stimulus or L-M cone stimulus (=1 for S cone)

% The LMS2RGB_Vpixx function will set it to the max possible contrast if 1.0 is beyond the gamut of the display.
cone_contrast = 1; % Which ever one we choose (S or L-M) will have this much contrast

% Set aside some of these parameters for saving:
ColParams.cone_contrast = cone_contrast;
ColParams.Scone = UseS_cone;

% If you want the 2 deg cone fundamentals, set to true:
Load2Deg = false; % otherwise, the 10 deg ones will be loaded.

% the number of trials PER EYE.
NumTrials = 3;

%If you're using the PROpixx or Viewpixx
UsingVP = true;
useHardwareStereo = false; %Setting to true works much better on the Viewpixx (+mac).

%Switch for whether you want the annulus superimposed over the dots:
DrawAnnulus = true;

% Set up unique file name to save the results in.
if UseS_cone
    ConeType = 'S_cone';
else
    ConeType = 'L-M_cone';
end
fName = fullfile('data', [SubjCode, '_isolum_', ...
    ConeType, '_', num2str(runN), '.mat']);

%Check that it doesn't already exist
if exist( fName, 'file' )
    userResponse = input( 'WARNING: File exists. Overwrite? Enter y or n: ', 's' );
    if ~strcmp( userResponse, 'y' )
        error('Aborting function' );
    end
end

%Choose the screen: it is usually the max screen no. available.
%Frustratingly, the Shuttle XPC (purchased June 2015) always seems to make the Vpixx display == 1. Not sure why, & can't seem to change it.
%So if we're on that machine, need to -1 from the screen number:
[~, CompName] = system('hostname'); %find out the computer name
if strncmpi(CompName, 'pspcawshuttle', length('pspcawshuttle')) ... %normal strcmp not working here, can't figure out why...
        && (length(Screen('Screens'))>1) %and there is more than 1 display connected...
    WhichScreen = max( Screen( 'Screens' ) )-1;
else
    WhichScreen = max( Screen( 'Screens' ) ); %should be right for any other machine!
end

screenRect = Screen('Rect',WhichScreen); %get the screen resolution.
centreX = screenRect(3)/2;
centreY = screenRect(4)/2;
RefreshRate = Screen('NominalFrameRate', WhichScreen);

%%%%-------------------------------%%%%
%         Load the spectra:           %
%%%%-------------------------------%%%%

% Obviously, we have different spectra for the Viewpixx and PROpixx. We load the appropriate ones
% depending on what machine we're running. If we're not on a Vpixx device (or not using it), we just load
% the Viewpixx spectra anyway.
% It is easier to load these here so they can be parsed to the LMS2RGB function.
% Here we also nominate the no. of pixels per degree, which also depends on the display.
if UsingVP
    Datapixx('Open');
    if Datapixx('IsViewpixx3D')  % If it's the Viewpixx3D
        load('C:\Users\Ryan Maloney\Documents\Calibration\Viewpixx_Processed_cal_data_2_4_2016.mat')
        %load('Viewpixx_Processed_cal_data_2_4_2016.mat')
        PPD = 37; %At 57cm viewing distance, there are 37 pixels/deg on the Viewpixx
        ColParams.Display = 'Viewpixx3D'; % Set aside display details
    elseif Datapixx('IsPropixx') % If it's the PROpixx projector
        % *** Here, we need to load the PROpixx spectra (don't need gamma because it's already linear):
        % NOTE: obviously, we would never load these when on the Viewpixx.
        %load('\\PSHome\Home\rm1380\My Documents\Calibration\PROpixx_Processed_cal_data_21_4_2016.mat')
        load('C:\Users\Ryan Maloney\Documents\Calibration\PROpixx_Processed_cal_data_21_4_2016.mat')
        PPD = 46.4; % Pixels per degree on the Propixx, also at 57 cm viewing distance
        ColParams.Display = 'PROpixx'; % Set aside display details
    end
else % If we're not using a Vpixx device, just load the Viewpixx spectra anyway (need spectra from somewhere)...
    load('C:\Users\Ryan Maloney\Documents\Calibration\Viewpixx_Processed_cal_data_2_4_2016.mat')
    %load('Viewpixx_Processed_cal_data_2_4_2016.mat')
    PPD = 37;
    ColParams.Display = 'NotVpixx'; % Set aside display details
end

% clear the java heap space if not on a mac:
if ~ismac
    jheapcl;
end

%%%%-------------------------------%%%%
%           Set up stimulus
%%%%-------------------------------%%%%

% Define the dot texture, a square-shaped sheet of dots.
%Make the texture the same size as the height of the screen
imsize = screenRect(4);

%Define dot density in dots/(deg^2):
dot_dens_per_deg2 = 1;

%compute number of dots in the total available area
num_dots = round(dot_dens_per_deg2 * (imsize/PPD)^2); %this many dots in the full dot field

%Just check whether num_dots is odd or even: (important for when contrast polarity is assigned)
%If odd, -1 to make it even
if mod(num_dots,2) %if odd, mod = 1
    num_dots = num_dots-1;
end

% specify dot size:
% Note, dot size needs to be larger for the cone-isolating stimuli than the black/white dots
dot_sigma_in_degrees = 0.2; %size of SD of the dot profile in degs/vis angle
dot_sigma = dot_sigma_in_degrees * PPD; %sigma in pixels
% Set the size of the circular dots drawn directly with 'DrawDots':
% there are 2.355 sigma under the curve of a Gaussian at FWHM.
% So the value assigned here approximates the width of our Gaussian dots.
circ_dot_size = dot_sigma * 2.355;

%NOTE: dotsize is simply half the length of the sides of a square patch of pixels that the dot profile is placed within.
%It is dot_sigma that really determines the size of the dots. Obviously, dotsize must be larger than dot_sigma.
%You also want it to be big enough that the convolution of the dots in the matrix allows for smooth transitions from dot colour to background colour
%for both black and white dots: ie so they are not 'aliased'
dotsize = round(dot_sigma * 10); %make the dots some multiple of sigma

%Define the minimum spacing between the PEAKS of the dot positions (in pixels):
SpatialJitter = round(0.5 * PPD); %+ dotsize/2;
%NOTE: add dot radius (ie dotsize/2) ensures the EDGES of the dots are separated by the minimum, but there may not be space in the matrix!

%determine random dot x and y positions, with minimum separation 'SpatialJitter'
%ie pseudo-random positioning to prevent dot overlap/clustering
ii=1
dot_pos = zeros(num_dots,2); %assign dot location matrix
tic
while ii <= num_dots
    
    if ii == 1 %set random coordinates for very first dot
        dot_pos(ii,:) = imsize.*rand(1,2);
        ii = ii+1
    else
        %set the next dot's random position
        dot_pos(ii,:) = imsize.*rand(1,2);
        %find the smallest distance (in pixels) between the current dot and any other existing dot
        idx = 1:ii-1; %index each existing dot except the current one
        d = min((dot_pos(idx,1)-dot_pos(ii,1)).^2 + (dot_pos(idx,2)-dot_pos(ii,2)).^2);
        d = sqrt(d);
        
        %Now if that distance is smaller than the minimum permitted distance, re-randomise the dot coordinates
        %This will continue until (at least) the minimum distance is met
        if d < SpatialJitter
            dot_pos(ii,:) = imsize.*rand(1,2);
        else %if that minimum distance is met, move on to the next dot
            ii = ii+1
        end
    end
end
toc

AdjustXBy = (screenRect(3) - screenRect(4))/2; %shift X dot positions by this amount to centre them, since image size <screenRect(3)
dot_pos(:,1) = dot_pos(:,1) + AdjustXBy;

%set up the raised cosine annular window.
%specify parameters for the annulus:
inrad = PPD * 1;% inner radius of annulus (in pixels), for fixation spot
outrad = PPD * 12/2; %outer radius of annulus (in pixels)
% define extent of spatial raised cosine at edge of aperture (in pixels)
% Here, we don't so much need the cosine ramping so just set it to something very small
% to give us a hard-edged annulus. Obviously this annulus doesn't need to be so complex
% if we're just using a circular edge but it's probably easier just to leave it this way.
cos_smooth = 2;
%This should plonk the window in the middle of the matrix, which is what we want
imsize2 = imsize*2; %double the texture size
x0 = (imsize2+1)/2;
y0 = (imsize2+1)/2;
J = ones(imsize2);
for (ii=1:imsize2)
    for (jj=1:imsize2)
        r2 = (ii-x0)^2 + (jj-y0)^2;
        if (r2 > outrad^2)
            J(ii,jj) = 0;
        elseif (r2 < inrad^2)
            J(ii,jj) = 0;
        elseif (r2 > (outrad - cos_smooth)^2)
            J(ii,jj) = cos(pi.*(sqrt(r2)-outrad+cos_smooth)/(2*cos_smooth))^2;
        elseif (r2 < (inrad + cos_smooth)^2)
            J(ii,jj) = cos(pi.*(sqrt(r2)-inrad-cos_smooth)/(2*cos_smooth))^2;
        end
    end
end

%%%%-------------------------------%%%%
%       Set up fixation
%%%%-------------------------------%%%%

% Set up the fixation cross or spot:
% This is drawn directly to the screen using Screen('FillRect')
% if you're using a cross instead:
crossWidth = 2;
crossHeight = 10;
fixationCross = fixation_cross(crossWidth,crossHeight,centreX,centreY);

% Make the fixation lock ring:
% We have an inner one around fixation and an outer one right on the edge of screen.
% These could probably be defined as a single texture (rather than 2) but I think that will complicate matters with the alpha-blending settings.
% (they are complicated enough already)
ringRadiusInner = PPD*0.5;                % ring surrounding fixation
ringRadiusOuter = screenRect(4)/2;        % outer edge (radius) of the ring: the edge of the screen
ringWidthInner = ringRadiusInner - PPD/4; % 1/4 of a degree thick
ringWidthOuter = ringRadiusOuter - PPD/3; % 1/3 of a degree thick

%Make the ring. It's in a 2*2 checkerboard pattern:
fixationRing = double(checkerboard(screenRect(4)/2,1) > 0.5);
%Define the ring:
xx = (1-imsize)/2:(imsize-1)/2;
[xx,yy] = meshgrid(xx,xx);
[~,r] = cart2pol(xx,yy);
% make the alpha mask for the rings, inner and outer.
ring_alphaInner = ((r>ringWidthInner+1) & (r<ringRadiusInner-1)); % Make the alpha mask a tiny bit thinner than the ring itself.
ring_alphaOuter = ((r>ringWidthOuter+1) & (r<ringRadiusOuter-1));

%%%%-------------------------------%%%%
%        Define response keys
%%%%-------------------------------%%%%

%Unify the keyboard names in case we run this on a mac:
KbName('UnifyKeyNames')

% Set up the keys for responses. If they seem a little awkward, it's because
% they are chosen to match the button response assignments on the Current Designs fORP-932 button
% response pad at the scanner, where the 4 buttons are registered as 1:4 on the keyboard (not numeric keypad)
% and triggers come through as '5'.
Increment = [KbName('4$'), KbName('4')]; % Increase theta
Decrement = [KbName('2@'), KbName('2')]; % Decrease theta
Accept = [KbName('1!'), KbName('1'), KbName('3#'), KbName('3')]; % Press these to accept current values.

% Finally, press a q at any time to quit the program.
RespQuit = KbName('q'); % q to quit.

try %Start a try/catch statement, in case something goes awry with the PTB functions
    
    %----------------------------
    % Set up the screen
    %----------------------------
    
    % initialization of the display
    AssertOpenGL;
    % Open PTB onscreen window: We request a 32 bit per colour component
    % floating point framebuffer if it supports alpha-blending. Otherwise
    % the system shall fall back to a 16 bit per colour component framebuffer:
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    % required for gamma correction through the PsychImaging pipeline:
    PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
    %Set the color range to be normalised between 0 and 1 (rather than 0-255):
    PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange', 1);
    
    % Open an on screen (grey) window and configure the imaging pipeline
    % Info about the 'blueline' mechanism for synching to the 3D glasses:
    % There seems to be a blueline generation bug on some OpenGL systems.
    % SetStereoBlueLineSyncParameters(windowPtr, windowRect(4)) corrects the
    % bug on some systems, but breaks on other systems.
    % We'll just disable automatic blueline, and manually draw our own bluelines!
    
    if useHardwareStereo
        [win, windowRect] = PsychImaging('OpenWindow', WhichScreen, 0.5, [], [], [], 1); %flag of 1 engages stereomode
        SetStereoBlueLineSyncParameters(win, windowRect(4)+10);
    else
        [win, windowRect] = PsychImaging('OpenWindow', WhichScreen, 0.5);
    end
    
    %Initialise the Vpixx device:
    if UsingVP        % Enable DATAPixx blueline support, and VIEWPixx scanning backlight for optimal 3D
        
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'UseDataPixx');
        Datapixx('Open');
        %The following commands are included in demos that apparently work for both the Viewpixx AND the PROpixx, though they seem specific to the Viewpixx...
        Datapixx('DisableVideoScanningBacklight');    % optionally, turn it off first, in case the refresh rate has changed since startup
        Datapixx('EnableVideoScanningBacklight');     % Only required if a VIEWPixx.
        Datapixx('EnableVideoStereoBlueline');
        Datapixx('SetVideoStereoVesaWaveform', 2);    % If driving NVIDIA glasses
        
        if Datapixx('IsViewpixx3D') %If it's the Viewpixx3D
            
            % Do the gamma correction. When using the Viewpixx/PROpixx through the PsychImaging pipeline, we
            % SHOULD NOT use Screen(‘LoadNormalizedGamma’) (see http://www.jennyreadresearch.com/research/lab-set-up/datapixx/)
            % The PROpixx device should have a linear lUT built in, but we will add this here for completeness.
            % The gamma values here were obtained following measurements (through the goggles) on
            % the Jaz Spectrometer taken 2/4/2016.
            % We will simply average the left and right eye's values, because they are so similar.
            
            gammaRed = mean([GammaValues(1,1,1), GammaValues(1,1,2)]);
            gammaGreen = mean([GammaValues(2,1,1), GammaValues(2,1,2)]);
            gammaBlue = mean([GammaValues(3,1,1), GammaValues(3,1,2)]);
            
            % We'll use the average of the right and left gamma values
            PsychColorCorrection('SetEncodingGamma', win, [1/gammaRed, 1/gammaGreen, 1/gammaBlue]);
            
            %Datapixx('EnableVideoLcd3D60Hz');
            Datapixx('DisableVideoLcd3D60Hz'); %=> According to Daniel B, disabling seems to give less crosstalk, bizarrely!
            subjectData.DisplayType = 'Viewpixx3D'; %set aside the device type for reference
            Datapixx('RegWr');
            
        elseif Datapixx('IsPropixx') %if it's the Propixx DLP projector
            
            subjectData.DisplayType = 'PROpixx'; %set aside the device type for reference
            Datapixx('SetPropixxDlpSequenceProgram',0); %set to normal RGB video processing for driving the LEDs & DLP MMDs
            %Datapixx('RegWr');
            
            %Modify the per-eye crosstalk on the PROpixx.
            %Apparently this cross-talk correction only works when using RB3D video mode,
            %where the red/blue channels contain the left/right eye greyscale images (which we are not using).
            %Datapixx('SetPropixx3DCrosstalkLR', 1);
            %Datapixx('SetPropixx3DCrosstalkRL', 1);
            Datapixx('RegWrRd'); %seem to need to do this after setting 'SetPropixxDlpSequenceProgram' to 0
        end
    end
    %No clue what the RegWr and RegWrRd commands are all about, but they are in the demos, so I've included them.
    
    %Define the 'blue line' parameters
    blueRectLeftOn   = [0,                 windowRect(4)-1, windowRect(3)/4,   windowRect(4)];
    blueRectLeftOff  = [windowRect(3)/4,   windowRect(4)-1, windowRect(3),     windowRect(4)];
    blueRectRightOn  = [0,                 windowRect(4)-1, windowRect(3)*3/4, windowRect(4)];
    blueRectRightOff = [windowRect(3)*3/4, windowRect(4)-1, windowRect(3),     windowRect(4)];
    
    HideCursor;
    %raise priority level:
    priorityLevel=MaxPriority(win); Priority(priorityLevel);
    %Query the screen refresh rate:
    ifi = Screen('GetFlipInterval',win); %in sec
    
    %Set the alpha-blending:
    %We want a linear superposition of the dots should they overlap:
    %Just like the Gabors in GarboriumDemo.m (see there for further info).
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE); %(Not sure this is actually working for the color dots)
    
    % We also want alpha-blending for smooth (anti-aliased) dots...
    %not sure how this will conflict with the above command
    %about the linear superposition of dots... but it doesn't seem to cause problems
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %==> definitely need this!!!
    
    %------------------------------------%
    %       Prepare LMS stimulus
    %------------------------------------%
    
    % Load the Stockman/Sharpe cone fundamentals:
    if Load2Deg % for 2 deg fundamentals
        load('C:\Users\Ryan Maloney\Documents\Colour\StockmanSharpe_2deg_cone_fundamentals_1nm.mat')
        %load('StockmanSharpe_2deg_cone_fundamentals_1nm.mat')
    else        % or for the 10 deg fundamentals
        load('C:\Users\Ryan Maloney\Documents\Colour\StockmanSharpe_10deg_cone_fundamentals_1nm.mat')
        %load('StockmanSharpe_10deg_cone_fundamentals_1nm.mat')
    end
    
    % Start off at a random point; some value close to the expected theta of pi/2.
    MaxIncrement = deg2rad(20); % Vary theta this much either side of pi/2
    ThetaIncr = deg2rad(1); % the amount we increment by each step
    % Initial theta:
    theta = (2*rand-1) * MaxIncrement + pi/2;
    
    if UseS_cone
        stimLMS.dir = [cos(theta)/sqrt(2), cos(theta)/sqrt(2), sin(theta)]; %the desired LMS cone activation
    else
        stimLMS.dir = [-sin(theta), cos(theta), 0]; %this eqn from Psykinematix; not sure if appropriate.
    end
    
    % Assign dot colour RGB values using LMS2RGB for both dot polarities for the first set of dots:
    % These are plugged directly into 'DrawDots' - no need to make textures.
    % We do it here mostly to initialise the variables, because they are reassigned again during presentation.
    % Note that because we are testing eyes separately, we need to present only grey to the
    % dot we're currently not testing, so that is simply RGB  = [0.5 0.5 0.5];
    % Do positive polarity dot:
    stimLMS.scale = cone_contrast; %the scale factor (akin to desired cone contrast, I think); max for s-cones is about 0.748 on the viewpixx
    [stimRGB, maxcontrast] = LMS2RGB_Vpixx (stimLMS, fundamentals, resampledSpectra);
    dotRGB(:,1) = stimRGB.dir * stimRGB.scale + 0.5; %scale & add the background to the obtained RGB triplet
    
    % Do negative polarity dot:
    stimLMS.scale = -cone_contrast;
    [stimRGB, maxcontrast] = LMS2RGB_Vpixx (stimLMS, fundamentals, resampledSpectra);
    dotRGB(:,2) = stimRGB.dir * stimRGB.scale + 0.5; %scale & add the background to the obtained RGB triplet
    
    % Set up separate variables for the 2 eyes:
    % We'll go in the order LRLRLR:
    dotRGB_L = dotRGB; % So start with the left eye.
    dotRGB_R = ones(3,2)*0.5; % Meaning with start with grey/blank for right eye.
    
    % Generate the annulus texture:
    AnnImages = 0.5.*ones(imsize2,imsize2,2); %specify RGB matrix of annulus, grey
    %specify Alpha channel of annulus, stays the same
    AnnImages(:,:,2) = J;
    annulus = Screen('MakeTexture',win,AnnImages,[],[],2);
    
    % Generate the Inner ring/fixation lock texture:
    ringMat(:,:,1) = fixationRing;
    ringMat(:,:,2) = ring_alphaInner;
    fixationRingTextureInner = Screen('MakeTexture',win,ringMat,[],[],2);
    
    % Generate the Outer ring/fixation lock texture:
    ringMat(:,:,1) = fixationRing;
    ringMat(:,:,2) = ring_alphaOuter;
    fixationRingTextureOuter = Screen('MakeTexture',win,ringMat,[],[],2);
    
    %%%%-------------------------------------------------%%%%
    %               Display the welcome screen
    %%%%-------------------------------------------------%%%%
    
    % Wait for the user to begin.
    Screen('TextFont',win, 'Arial');
    Screen('TextSize',win, 24);
    if useHardwareStereo
        Screen('SelectStereoDrawBuffer', win, 0);   % flag of 0= left eye
        DrawFormattedText(win,['Welcome ' SubjCode , ...
            '. \n \nPress '' Red '' or '' 4 '' to increase luminance,', ...
            '\n \nor press '' Yellow '' or '' 2 '' to decrease.', ...
            '\n \nPress '' Green '' or '' Blue '' (1 or 3) to accept the settings and move on.', ...
            '\n \nPress '' q '' on the keyboard to quit at any time.', ...
            '\n \nPress any button/key to begin.'], ...
            'center', 'center', 0);
        
        Screen('SelectStereoDrawBuffer', win, 1);   % flag of 1= right eye
        DrawFormattedText(win,['Welcome ' SubjCode , ...
            '. \n \nPress '' Red '' or '' 4 '' to increase luminance,', ...
            '\n \nor press '' Yellow '' or '' 2 '' to decrease.', ...
            '\n \nPress '' Green '' or '' Blue '' (1 or 3) to accept the settings and move on.', ...
            '\n \nPress '' q '' on the keyboard to quit at any time.', ...
            '\n \nPress any button/key to begin.'], ...
            'center', 'center', 0);
        Screen('Flip', win); % , [], [], [], 1);
    else
        DrawFormattedText(win,['Welcome ' SubjCode , ...
            '. \n \nPress '' Red '' or '' 4 '' to increase luminance,', ...
            '\n \nor press '' Yellow '' or '' 2 '' to decrease.', ...
            '\n \nPress '' Green '' or '' Blue '' (1 or 3) to accept the settings and move on.', ...
            '\n \nPress '' q '' on the keyboard to quit at any time.', ...
            '\n \nPress any button/key to begin.'], ...
            'center', 'center', 0);
        Screen('Flip', win); %, [], [], [], 1);
    end
    
    WaitSecs(0.2)
    KbCheck(); % take a quick KbCheck to load it now & flush any stored events
    
    % Wait for user response to continue...
    ButtonPressed = 0;
    while ~ButtonPressed
        [KeyIsDown, ~, keyCode] = KbCheck();
        if KeyIsDown
            if keyCode(RespQuit)
                ForcedQuit = true
                ExitGracefully(UsingVP, ForcedQuit)
                %if any other button on the keyboard has been pressed
            else
                ButtonPressed = 1;
            end
        end
    end
    
    %%%%---------------------------------------------------------------%%%%
    %               Set up dot colour & location indices
    %%%%---------------------------------------------------------------%%%%
    
    % We need to set up contrasts & timings (dot lifetimes) for each trial and each frame.
    % We will assume we're running on the Vpixx device with a 120 Hz frame rate (so 60 Hz per-eye).
    % PER-EYE frames 1:fperiod/2 are one colour, per-eye frames fperiod/2:fperiod are the other.
    % Each dot will have a fperiod/2 per-eye frame lifetime, so we must index the rows of RGB appropriately
    % we can do this using fperiod = n (ie n frame cycle) and, [mod(f,fperiod) < (fperiod/2)] +1
    % Note that the 'blank' periods for each alternate eye are included in the total dot lifetime:
    % ie we compute dot lifetime on the per-eye frame rate
    fperiod = 8; % the PER EYE frames per colour CYCLE (ie YVYVYV ...)
    % Compute & store the dot lifetime.
    ColParams.DotLifetimeMS = 1000/(RefreshRate/2) * fperiod/2; %=> the dot lifetime for a single colour.
    
    f = 0; % this value increases with each iteration, but on a PER-EYE basis only
    missedFrames = 0;
    vbl = Screen('Flip',win); % sync vbl to start time
    ForcedQuit = false;
    
    %%%%------------------------------------------------%%%%
    %               Run the experiment!
    %%%%------------------------------------------------%%%%
    
    for trial = 1:NumTrials * 2 % run for the number of trials * 2 eyes
        
        trial % print out current trial
        IsoAccept = 0; % to control acceptance of the current settings/move to next trial
        MoveToNext = 0;
        
        %Loop across stimulus frames:
        while ~IsoAccept
            
            % Select left-eye image buffer for drawing:
            if useHardwareStereo
                Screen('SelectStereoDrawBuffer', win, 0);
            end
            
            %%%%------------------------------------------------%%%%
            %               Draw left eye stimulus:
            %%%%------------------------------------------------%%%%
            
            %Draw dots:
            Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE); %turn on alpha blending for lin superposition: see GarboriumDemo.m
            Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('DrawDots', win, dot_pos', ...                 % the dot xy coords
                circ_dot_size, ...                                % the dot width
                dotRGB_L(:,(mod(f,fperiod) < (fperiod/2))+1), ... % dot colour, for current frame
                [],2);                                            % Final flag of 2 gives us anti-aliased dots.
            %Screen('DrawTextures',win,DotsIdx(mod(f,6)+1,:),[],dstRects',[],[],1)
            
            %Superimpose the annulus:
            if DrawAnnulus
                Screen('Blendfunction', win, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                Screen('DrawTexture',win,annulus);
                Screen('Blendfunction', win, GL_ONE, GL_ZERO);
            end
            
            %Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
            Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
            Screen('DrawTexture',win, fixationRingTextureInner);
            Screen('DrawTexture',win, fixationRingTextureOuter);
            %Draw the black fixation cross:
            Screen('FillRect',win,[0 0 0],fixationCross);
            Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); %Done drawing so flip alpha back again
            Screen('BlendFunction', win, GL_ONE, GL_ZERO);
            %Draw blue lines:
            Screen('FillRect', win, [0, 0, 1], blueRectLeftOn);
            Screen('FillRect', win, [0, 0, 0], blueRectLeftOff);
            
            % Select right-eye image buffer for drawing:
            if useHardwareStereo
                Screen('SelectStereoDrawBuffer', win, 1);
            else %Not sure what would happen if this 'else' was actually true. I suspect something would go wrong, as stim would be presented twice.
                %But this is more or less how the Vpixx people do it in their 'DatapixxImagingStereoDemo'
                Screen('DrawingFinished', win);
                [vbl , ~, ~, missed] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); %update display on next refresh (& provide deadline)
                %f = f+1
                if missed > 0
                    missedFrames = missedFrames + 1;
                end
            end
            
            %%%%------------------------------------------------%%%%
            %               Draw right eye stimulus:
            %%%%------------------------------------------------%%%%
            
            %Draw dots:
            Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE); %turn on alpha blending for lin superposition: see GarboriumDemo.m
            Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('DrawDots', win, dot_pos', ...                 % the dot xy coords
                circ_dot_size, ...                                % the dot width
                dotRGB_R(:,(mod(f,fperiod) < (fperiod/2))+1), ... % dot colour, for current frame
                [],2);
            %Screen('DrawTextures',win,DotsIdx(mod(f,6)+1,:),[],dstRects',[],[],1)
            
            %Superimpose the annulus:
            if DrawAnnulus
                Screen('Blendfunction', win, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                Screen('DrawTexture',win,annulus);
                Screen('Blendfunction', win, GL_ONE, GL_ZERO);
            end
            
            %Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
            Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
            Screen('DrawTexture',win, fixationRingTextureInner);
            Screen('DrawTexture',win, fixationRingTextureOuter);
            %Draw the fixation cross:
            Screen('FillRect',win,[0 0 0],fixationCross);
            Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); %Done drawing so flip alpha back again
            Screen('BlendFunction', win, GL_ONE, GL_ZERO);
            %Draw blue lines:
            Screen('FillRect', win, [0, 0, 1], blueRectRightOn);
            Screen('FillRect', win, [0, 0, 0], blueRectRightOff);
            
            Screen('DrawingFinished', win);
            
            [vbl, ~, ~, missed] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); %update display on next refresh (& provide deadline)
            
            f = f+1; %increment counter for next PER-EYE frame
            %keep record of any missed frames:
            if missed > 0
                missedFrames = missedFrames + 1;
            end
            
            %%%%------------------------------------------------%%%%
            %               Check for response
            %%%%------------------------------------------------%%%%
            
            % Check for keyboard response
            [keyIsDown, ~, keyCode] = KbCheck();
            ValidResp = 0;
            if keyIsDown %If they've pressed a button on keyboard,
                WaitSecs(0.05);
                % determine if it is valid and proceed
                
                %If they respond 'Increase'
                if any(keyCode(Increment))
                    theta = theta + ThetaIncr
                    IsoAccept = 0;
                    % If they reach the outer limits of possible theta.
                    % Note that they will get no where once they've reached this point
                    % if they continue in the same direction.
                    if theta > (MaxIncrement+pi/2)
                        theta = MaxIncrement + pi/2;
                        
                        % Play audio feedback to indicate upper end of the range:
                        Beeper(2000,10, 0.05)
                        Beeper(2000,10, 0.05)
                    end
                    ValidResp = 1;
                    
                    % If they respond 'Decrease'
                elseif any(keyCode(Decrement))
                    IsoAccept = 0;
                    theta = theta - ThetaIncr
                    
                    % As above, set the limit on how far below pi/2 they can go.
                    if theta < (pi/2 - MaxIncrement)
                        theta = pi/2 - MaxIncrement;
                        Beeper(500,10, 0.05)
                        Beeper(500,10, 0.05)
                    end
                    ValidResp = 1;
                    
                    % If they accept the current settings:
                elseif any(keyCode(Accept))
                    
                    % Set aside the chosen theta value:
                    Col_results.theta_vals(trial) = theta;
                    
                    %Reset theta to a new start value for next trial:
                    theta = (2*rand-1) * MaxIncrement + pi/2;
                    f = 0; %reset f to go back to 1st frame.
                    ValidResp = 1;
                    % MoveToNext will break the current while loop & move to next trial, but only after the colour
                    % values have been randomly re-assigned. It also controls which eye is currently being tested.
                    MoveToNext = 1;
                    
                elseif keyCode(RespQuit)
                    % break out of program if 'q' is pressed
                    ForcedQuit = true;
                    ExitGracefully(UsingVP, ForcedQuit)
                    
                end
                
                % Only re-evaluate the LMS settings if a valid response was made.
                % This includes an increment or decrement, or an 'accept'.
                % But even if it's been accepted we must still prepare the stimulus for the next trial.
                % A random value for theta has already been assigned above.
                if ValidResp
                    % Re-evaluate the stimulus:
                    tic
                    if UseS_cone
                        stimLMS.dir = [cos(theta)/sqrt(2), cos(theta)/sqrt(2), sin(theta)]; %the desired LMS cone activation
                    else
                        stimLMS.dir = [-sin(theta), cos(theta), 0]; %this eqn from Psykinematix; not sure if appropriate.
                    end
                    stimLMS.scale = cone_contrast; %the scale factor (akin to desired cone contrast, I think); max for s-cones is about 0.748 on the viewpixx
                    [stimRGB, maxcontrast] = LMS2RGB_Vpixx (stimLMS, fundamentals, resampledSpectra);
                    dotRGB(:,1) = stimRGB.dir * stimRGB.scale + 0.5; %scale & add the background to the obtained RGB triplet
                    
                    % Do negative polarity dot:
                    stimLMS.scale = -cone_contrast;
                    [stimRGB, maxcontrast] = LMS2RGB_Vpixx (stimLMS, fundamentals, resampledSpectra);
                    dotRGB(:,2) = stimRGB.dir * stimRGB.scale + 0.5 %scale & add the background to the obtained RGB triplet
                    toc
                    
                    % Update RGB for appropriate eye. Note this only occurs if a valid response has been made,
                    % so they are not updated on every frame, which would be silly because there's no need
                    % to update them. If we move to the next trial they are updated again but that doesn't matter.
                    
                    
                    % If the following == 1, it means we're on the R eye, so RGB goes into R, grey into L.
                    % otherwise, if 0, we do the opposite
                    if mod(trial+1,2) % if yes, then right eye dots
                        dotRGB_L = ones(3,2)*0.5;
                        dotRGB_R = dotRGB;
                    else              % if no, then left eye dots
                        dotRGB_L = dotRGB;
                        dotRGB_R = ones(3,2)*0.5;
                    end
                    
                    % If the current settings have been accepted, prepare for next trial & move on
                    if MoveToNext
                        
                        % First, set up a new set of dot positions for the next trial, just
                        % as we did above:
                        ii=1;
                        dot_pos = zeros(num_dots,2); %assign dot location matrix
                        while ii <= num_dots
                            if ii == 1 %set random coordinates for very first dot
                                dot_pos(ii,:) = imsize.*rand(1,2);
                                ii = ii+1;
                            else
                                %set the next dot's random position
                                dot_pos(ii,:) = imsize.*rand(1,2);
                                %find the smallest distance (in pixels) between the current dot and any other existing dot
                                idx = 1:ii-1; %index each existing dot except the current one
                                d = min((dot_pos(idx,1)-dot_pos(ii,1)).^2 + (dot_pos(idx,2)-dot_pos(ii,2)).^2);
                                d = sqrt(d);
                                %Now if that distance is smaller than the minimum permitted distance, re-randomise the dot coordinates
                                %This will continue until (at least) the minimum distance is met
                                if d < SpatialJitter
                                    dot_pos(ii,:) = imsize.*rand(1,2);
                                else %if that minimum distance is met, move on to the next dot
                                    ii = ii+1;
                                end
                            end
                        end
                        %shift X dot positions by this amount to centre them, since image size <screenRect(3)
                        AdjustXBy = (screenRect(3) - screenRect(4))/2;
                        dot_pos(:,1) = dot_pos(:,1) + AdjustXBy;
                        
                        % if the following == 1, it means we're currently on L eye, about to move to R eye,
                        % So assign the RGB values accordingly. 0.5 (blank/grey) for non-tested eye.
                        % Note that these are the opposite conditons to the above scenario if still within the same trial.
                        if mod(trial,2)
                            dotRGB_L = ones(3,2)*0.5;
                            dotRGB_R = dotRGB;
                        else
                            dotRGB_L = dotRGB;
                            dotRGB_R = ones(3,2)*0.5;
                        end
                        
                        IsoAccept = 1 % Break the current while loop, move to next
                        
                    end
                end     % end of actions following valid response (LMS updating)
            end         % end of checking for response
        end             % end of loop across stimulus frames
        
        % Blank the screen briefly between trials:
        Screen('Flip', win);
        WaitSecs(1.5)
    end                 % end of loop across trials
    
    % Compute the mean theta for the 2 eyes:
    Col_results.MeanThetas(1) = mean(Col_results.theta_vals(1:2:end)); % mean theta: left eye
    Col_results.MeanThetas(2) = mean(Col_results.theta_vals(2:2:end)); % mean theta: right eye
    
    % Helper: compute frequency of flicker.
    % This is the frequency by which the dots flick from one colour to the other.
    % Note that dot duration includes the 'blank' periods while the other eye is stimulated.
    ColParams.FlickerFrequency = 1000/(ColParams.DotLifetimeMS * 2);
    
    % Also include the max possible contrast given the display gamut:
    ColParams.MaxPossibleContrast = maxcontrast;
    
    % Save the results:
    save(fName, 'Col_results', 'ColParams')
    
    %Print out the results:
    fprintf('\nThetas  = %.3f', Col_results.theta_vals);
    fprintf('\n \nMean theta, left eye = %.3f \n \nMean theta, right eye = %.3f\n', Col_results.MeanThetas(1), Col_results.MeanThetas(2));
    fprintf('\nMissed frames = %.3f \n', missedFrames)
    
    % Close everything down & exit:
    ExitGracefully (UsingVP, ForcedQuit)
    
catch MException
    
    rethrow (MException)
    psychrethrow(psychlasterror)
    ExitGracefully (UsingVP, ForcedQuit)
    error('Error!')
    
end %End of try/catch statement
end % End of main function.

% The end, Sub-functions follow...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rect = fixation_cross(width,height,centrex,centrey)
% $Id: fixation_cross.m 4 2007-06-23 11:13:19Z damienm $
%Make a small fixation cross to display on screen.
%Width and height are in pixels.
%centrex & centrey give the x & y screen coordinates where you want the cross centred.

rect = zeros(4,2);

width = width/2;
height = height/2;

rect(1,1) = -height;
rect(2,1) = width;
rect(3,1) = height;
rect(4,1) = -width;

rect(1,2) = -width;
rect(2,2) = height;
rect(3,2) = width;
rect(4,2) = -height;


rect(1:2:4,:) = rect(1:2:4,:) + centrex;
rect(2:2:4,:) = rect(2:2:4,:) + centrey;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExitGracefully (UsingVP, ForcedQuit)
%...need to shut everything down here...

% turn off the prioritisation:
Priority( 0 ); % restore priority

if UsingVP        % close down the ViewPixx or ProPixx
    Datapixx('DisableVideoScanningBacklight');
    if Datapixx('IsViewpixx3D')
        Datapixx('DisableVideoLcd3D60Hz');
    end
    Datapixx('RegWr');
    %Datapixx('Close'); %closing it here might cause it to crash?
end

% Close down the screen:
Screen('CloseAll')

Datapixx('Close'); % closing the Datapixx here (after closing the screen) might stop it from crashing

% Bring back the mouse cursor:
ShowCursor();

% announce to cmd window if the program was aborted by the user
if ForcedQuit
    error('You quit the program!')
end

end



