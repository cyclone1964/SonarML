%generateSamples - make a range doppler map with a target it in
%
% generateSamples(FileName,{PropertyList}) will generate a range
% doppler map with a target in it. It does this by making a time
% series with background noise, a reverb ridge, and a single
% pulse. The properties are:
%
% Band - band of sonar (2 x 1 Hz) [18000 22000]
% Baffling - amount of backplane baffling  (0 = none, nepers)
% Aperture - the aperture of the array in wavelengths
% Steering - the steering direction of the array
% FrameSize - size of analysis FFT (dimensionless) (256)
% Environment - environment ot use
% TargetRanges - Range vector of targets (NumTargetsx3, m) (500)
% SourceLevel - source level of emitter (dB) (200)
% CycleLength - length of cycle (s) (2)
% PlatformSpeed - speed of platform (m/s) (0)
% TargetStrength - target strength (dB) (20)
% TargetVelocity - doppler of target (m/s) (5)
% BackgroundLevel - background noise level (dB/rootHz) (100)
% VolumeReverbLevel - scattering level of volumereverb (dB) (-60)
% BoundaryReverbAdjustment - dB adjustment for Surface (1) and Bottom (2) (0)
%
% The model assumes an omnidirectional transmitter and a platform
% moving in the +X direction. This is all done at baseband with no
% absorption. 
% 
function [Samples,Properties] = generateSamples(varargin)

% Set default properties and then modify for input list
Properties.Band = [18000 20000];
Properties.Baffling = [];
Properties.Aperture = [];
Properties.Steering = [1 0 0]';
Properties.FrameSize = [];
Properties.Environment = initializeEnvironment;
Properties.SourceLevel = 220;
Properties.PulseLength = 0.50;
Properties.CycleLength = 1;
Properties.PlatformSpeed = 0;
Properties.PlatformDepth = 300;
Properties.TargetStrengths = [];
Properties.TargetPositions = [];
Properties.TargetVelocities = [];
Properties.BackgroundLevel = 40;
Properties.VolumeReverbAdjustment = 0;
Properties.BoundaryReverbAdjustment = [0 0]';

Properties = setProperties(Properties,varargin{:});

% Compute the sample rate, the bin width, the number of
% Samples, the number of frames, etc
SampleRate = diff(Properties.Band);
SoundSpeed = Properties.Environment.WaterColumn.SoundSpeed;

% If no frame size is given, let's estimate it from the length of the pulse
% and the sample rate;
if (isempty(Properties.FrameSize))
  Properties.FrameSize = ...
    2^(nextpow2(SampleRate * Properties.PulseLength)-1);
end
NumBins = Properties.FrameSize;

% Round length up to integer number of frames
HzPerBin = SampleRate/NumBins;
SecPerFrame = NumBins/SampleRate;
CenterFreq = mean(Properties.Band);
SamplesPerFrame = NumBins;
NumFrames = ceil(Properties.CycleLength/SecPerFrame);
NumSamples = NumFrames * NumBins;

% Extract a few things from the envioronment
SoundSpeed = mean(Properties.Environment.WaterColumn.SoundSpeed);
WaterDepth = Properties.Environment.Bottom.Depth;

% Initialize with the noise
NoiseLevel = Properties.BackgroundLevel - 10 * log10(SampleRate);
Noise = 10^(NoiseLevel/20) * ...
        (randn(NumSamples,1) + i * randn(NumSamples,1));

% Now add the signals after computing the doppler
NumTargets = max([size(Properties.TargetPositions,2) 
                  size(Properties.TargetVelocities,2)
                  length(Properties.TargetStrengths)]);

Signal = zeros(size(Noise));
for TargetIndex = 1:NumTargets
  
  Range = norm(Properties.TargetPositions(:,TargetIndex));
  Direction = Properties.TargetPositions(:,TargetIndex)/Range;
  SourceDoppler = Properties.PlatformSpeed * ([1 0 0] * Direction);
  TargetDoppler = Properties.TargetVelocities(:,TargetIndex)' * Direction;
  TotalDoppler = ...
      (1+SourceDoppler/SoundSpeed) * ...
      (1+TargetDoppler/SoundSpeed);

  PulseLevel = Properties.SourceLevel - ...
      40 * log10(Range) + ...
      Properties.TargetStrengths(TargetIndex);

  Beam = computeBeamResponse(Direction'*Properties.Steering,...
                             Properties.Aperture) .* ...
         computeBaffling(Direction(1),Properties.Baffling);

  PulseLevel = PulseLevel + 20*log10(Beam);

  NumPulseSamples = ...
    floor(Properties.PulseLength/TotalDoppler * SampleRate);
  BasebandFrequency = 2*pi*(TotalDoppler-1) *  CenterFreq/SampleRate;
  Phase = (1:NumPulseSamples)*BasebandFrequency;
  Pulse = (cos(Phase') + i * sin(Phase')) * 10^(PulseLevel/20);

  % Window that baby
  Pulse = Pulse .* blackman(length(Pulse));

  % Add to the signal
  Indices = round(2 * Range * SampleRate/ SoundSpeed) + (1:NumPulseSamples);
  Indices = Indices(Indices < NumSamples);
  Signal(Indices) = Signal(Indices) + Pulse(1:length(Indices));
end

% Now let's make some volume reverb. To do this we need to
% integrate the beam pattern into a spectrum.
AngleIncrement = pi/180;
Bearings = -(pi-AngleIncrement/2):AngleIncrement:(pi-AngleIncrement/2); 
Elevations = ...
    -(pi/2-AngleIncrement/2):AngleIncrement:(pi/2-AngleIncrement/2); 

% Now form integration points for all of them
[AllBearings,AllElevations] = meshgrid(Bearings,Elevations);

% Make a beam pattern
Directions = ...
    computeDirection([zeros(size(AllBearings(:))) ...
                    AllElevations(:) ...
                    AllBearings(:)]');

Beam = computeBeamResponse(Properties.Steering'*Directions, ...
                           Properties.Aperture) .* ...
         computeBaffling(Directions(1,:),Properties.Baffling);       

% And the differential area
DifferentialArea = cos(AllElevations(:)) * AngleIncrement^2;

% Now compute the frequency bin for all of them.
Temp = cos(AllElevations(:)) .* cos(AllBearings(:));
Bins = round(CenterFreq * ...
             (2*Properties.PlatformSpeed*Temp/SoundSpeed)/ ...
             HzPerBin) + NumBins/2;

% And integrate them into a spectrum
Spectrum = accumarray(Bins,DifferentialArea .* Beam(:));
if (length(Spectrum) < NumBins) 
  Spectrum(NumBins) = 0;
end
Spectrum = sqrt(fftshift(Spectrum));

% This is a hack until I can get the real window
Phase = pi/2 * (0:(SamplesPerFrame-1))/(SamplesPerFrame-1);
Window = generateMMWindow(NumBins);

% OK, having now done that, we make a series of frames of complex
% gaussian noise shaped by that
VolumeReverb = zeros(NumSamples,1);
LastFrame = [];
SampleIndex = 1;
for FrameIndex = 1:NumFrames
  Frame = NumBins * sqrt(2) * ...
    ifft(Spectrum .* (randn(NumBins,1) + i * randn(NumBins, 1)));
  Indices = (1:SamplesPerFrame) + SampleIndex-1;
  VolumeReverb(Indices) = Frame .* Window.Front;
  if (~isempty(LastFrame))
    VolumeReverb(Indices) = VolumeReverb(Indices) + LastFrame .* Window.Back;
  end
  LastFrame = Frame;
  SampleIndex = SampleIndex + NumBins;
end
Indices = 1:SamplesPerFrame;
VolumeReverb(Indices) = VolumeReverb(Indices) + LastFrame .* Window.Back;

% Now to first approximation let's just scale that by the two way
% prop loss
Times = ((1:length(VolumeReverb))-1)'/SampleRate + Properties.PulseLength;
Ranges = Times * SoundSpeed/2;
PropagationLoss = 1 ./ Ranges;
Level = Properties.SourceLevel + ...
        Properties.Environment.WaterColumn.ScatteringStrength + ...
        Properties.VolumeReverbAdjustment + ...
        10 * log10(HzPerBin) + ...
        10 * log10(SoundSpeed * Properties.PulseLength/2);
Scale = 10^(Level/20) * PropagationLoss;
VolumeReverb = VolumeReverb .* Scale;

% Now we have to make boundary reverb for the two boundaries. This
% assumes a constant scattering strength, which is of course a lie,
% but we don't have the curves around so that is that for now. 

% First, get times for the two boundaries where reverb actually hit
FrameTimes = ((1:NumFrames)-1)*SecPerFrame + ...
    Properties.PulseLength;

BoundaryDistances = [Properties.PlatformDepth
                     WaterDepth - Properties.PlatformDepth];

% Angles for the integration (this one in cylindrical coordinates
Angles = (1:360)-0.5;
BoundaryReverb = zeros(NumSamples,2);

for BoundaryType = 1:2

  % Find those frames for each of the two boundaries that contribute to
  % the signal
  Frames = ...
      find(FrameTimes > ...
           2*BoundaryDistances(BoundaryType)/SoundSpeed);
  
  % And now the Times of thoes frames
  FrontTimes = FrameTimes(Frames);
  BackTimes = FrontTimes + Properties.PulseLength;
  
  % These are the HORIZONTAL ranges, needed to form the scattering annuli
  FrontRanges = ...
      sqrt((FrontTimes*SoundSpeed/2).^2 - ...
           BoundaryDistances(BoundaryType)^2);
  BackRanges = ...
      sqrt((BackTimes*SoundSpeed/2).^2 - ...
           BoundaryDistances(BoundaryType)^2);
  FrameRanges = BackRanges - FrontRanges;
  
  % And the elevations to each frame
  Elevations = atan2d(BoundaryDistances(BoundaryType),FrontRanges);
  
  % The last frame for the MM windowing
  LastFrame = [];
  
  % Now do all the frames
  for FrameIndex = 1:(length(Frames)-1)
    
    % Compute the directions to the scatterin point
    Directions = ...
        computeDirection([zeros(size(Angles))
                        Elevations(FrameIndex) * ones(size(Angles))*pi/180
                        Angles*pi/180]);
    
    % Compute the frequency bins
    Bins = round(CenterFreq * ...
                 (2*Properties.PlatformSpeed/SoundSpeed) *  ...
                 [1 0 0] * Directions / HzPerBin) + NumBins/2;
    
    % Now compute the propagation loss, this time two way since we
    % don't gain the entire thing
    PropagationLoss = -40*log10(FrontTimes(FrameIndex)*SoundSpeed/2);
    Angle = abs(Elevations(FrameIndex)*pi/180);
    switch(BoundaryType)
      case 1
        Scattering = ...
            computeSurfaceScatteringStrength(CenterFreq, ...
                                             Angle, ...
                                             Properties.Environment);
      case 2
        Scattering = ...
            computeBottomScatteringStrength(CenterFreq, ...
                                            Angle, ...
                                            Properties.Environment);
    end

    Level = Scattering + ...
      Properties.SourceLevel + ...
      PropagationLoss + ...
      Properties.BoundaryReverbAdjustment(BoundaryType) + ...
      10*log10(HzPerBin) + ...
      10*log10(FrontRanges(FrameIndex)) + ...
      10*log10(FrameRanges(FrameIndex));
    
    Energy = 10.^(Level/10) * min(diff(Angles)*pi/180);
    
    % Now if there is an aperture supplied, we need to apply it
    Energy = Energy * ...
             computeBeamResponse(Properties.Steering'*Directions, ...
                                 Properties.Aperture) .*  ...
             computeBaffling(Directions(1,:),Properties.Baffling);                     

    % Now integrate into a spectrum
    Spectrum = accumarray(Bins',Energy);
    if (length(Spectrum) < NumBins)
      Spectrum(NumBins) = 0;
    end
    Spectrum = sqrt(fftshift(Spectrum));
    
    % Now make the frame
    Frame = NumBins * sqrt(2) * ...
      ifft(Spectrum .* (randn(NumBins,1) + i * randn(NumBins,1)));
    
    % And place it into the reverberation
    Indices = (1:SamplesPerFrame) + ((Frames(FrameIndex)-1)*SamplesPerFrame);
    if (isempty(LastFrame)) 
      BoundaryReverb(Indices,BoundaryType) = Window.Front .* Frame;
    else
      BoundaryReverb(Indices,BoundaryType) = ...
          Window.Back .* LastFrame + Window.Front .* Frame;
    end
    LastFrame = Frame;
  end
end

% Sum it all up
Samples = Noise + Signal+ ...
          VolumeReverb + BoundaryReverb(:,1) + BoundaryReverb(:,2);
      
%computeBeamResponse - compute a nominal beam pattern
%
% Beams = computeBeamPattern(DirCos) computes a beam pattern for a
% nominal cylindrical piston with a backplane baffle
% where the input is the cosine of the off-axis
% angle. 
function Beams = computeBeamResponse(DirCos, Aperture, Baffling)

% If there is no aperture, it is omnidirectional
if (isempty(Aperture) || Aperture <= 0)
    Beams = ones(size(DirCos));
    return;
end

% Do the circular piston part first, clipping it when it hits endfire
Args = 2*pi*Aperture*max(0,sqrt(1-DirCos.^2));
Args = max(1e-10,Args);
Beams = 2*besselj(1,Args)./Args;
Beams = Beams .* Beams;

%computeBaffling - compute the baffling of the array
%
% Baffling= computeBaffling(DirCos,Baffling) will return the
% backplane baffling of the array
%
function Baffling = computeBaffling(DirCos,Baffling)

% Now attempt to put a backplane baffling, starting 
% 0.1 in front of the plane
if (isempty(Baffling))
    Baffling = ones(size(DirCos));
    return;
end

Argument = min(0,DirCos-0.1);
Baffling = exp(Baffling * Argument);

