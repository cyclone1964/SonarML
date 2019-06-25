% This script makes a much harder data set, with the target 25 dB
% quieter than the simple set
Plot = true;

% THe first thing to do is reset the random number generator
rng('default');

% This is the number of frequency bins in a frame
NumBins = 32;

% For this simulation, the band is between 18 and 22 k
Band = [18000 22000];

% That means that the complex baseband sample rate is the bandwidth
SampleRate = diff(Band);
CenterFrequency = mean(Band);
PulseLength = NumBins/SampleRate;
SecondsPerFrame = NumBins/SampleRate;

% We want cycles that are 1 second long, so it's this many samples,
% but we round up to an integer number of frames
SamplesPerCycle = SampleRate;
FramesPerCycle = ceil(SamplesPerCycle/NumBins);
SamplesPerCycle = FramesPerCycle*NumBins;
MetersPerFrame = 750*SecondsPerFrame;

% Now, the number of Cyles. We want a total of 8k training
% spectra spread across these cycles.
NumCycles = ceil((8*1024)/FramesPerCycle);

% In order for the training to work, we want multiple targets per
% cycle. Otherwise there is a large imbalance in the
% target/notarget classes. however, because each target can exist
% in two frames, we need to reduce the number from /2 to /4 so they
% don't overlap
TargetsPerCycle = floor(FramesPerCycle/4);

% Now, let's make parameters for each of those cycles. We have the
% sonar speed uniform between -3 and 3, and the Dopplers uniform
% between -10 and 20
Speeds = 12*rand(TargetsPerCycle,NumCycles)-6;
Dopplers = 20 * rand(TargetsPerCycle,NumCycles)-10;

% Now, we need the range fixed to specific frames, so let's do that
% now, keeping them out of the first two and the last two frames
RangeBins = zeros(size(Dopplers));
for CycleIndex =1:NumCycles
    RangeBins(:,CycleIndex) = randperm(FramesPerCycle-4, TargetsPerCycle)+2;
end
Ranges = MetersPerFrame * (RangeBins-1);

% Open the file and write the header
FID = fopen('QuietData.csv','w');

% This environment is meant to be pretty clean
SoundSpeed = 1500;

WaterColumn = initializeWaterColumn('SoundSpeed',SoundSpeed, ...
                                    'ScatteringStrength',-75);
SurfaceBoundary = ...
    initializeSurfaceBoundary('SeaState',1,'WaveHeight',1,'WindSpeed',3);
BottomBoundary = ...
    initializeBottomBoundary('Depth',150,'Porosity',0.1,'GrainSize',1);
Environment = initializeEnvironment('SurfaceBoundary',SurfaceBoundary, ...
                                    'BottomBoundary',BottomBoundary, ...
                                    'WaterColumn',WaterColumn);

% Now lets go through every cycle and genrate the signals, make the
% range doppler map, and write it out. Let's also save the images
% of the cycles.
for CycleIndex = 1:length(Ranges)

    Speed = Speeds(CycleIndex);
    Doppler = Dopplers(CycleIndex);
    
    Angles = 2*pi*randn(1,TargetsPerCycle);
    Positions = [cos(Angles);sin(Angles);zeros(size(Angles))];
    Positions = repmat(Ranges(:,CycleIndex)',3,1) .* Positions;
    
    [Samples,Properties] = ...
        generateSamples('Band',Band, ...
                        'PulseLength',PulseLength, ...
                        'TargetStrengths',20*ones(size(Positions)), ...
                        'TargetPositions',Positions, ...
                        'TargetDopplers',Dopplers(:,CycleIndex), ...
                        'PlatformSpeed',Speed, ...
                        'PlatformDepth',400, ...
                        'Environment',Environment, ...
                        'VolumeReverbAdjustement',-20, ...
                        'BoundaryReverbAdjustment',[-100 -100]);
    
    % Make a spectrogram
    SoundSpeed = SoundSpeed;

    [S,Frequencies,Times,P] = spectrogram(Samples,NumBins,0,NumBins,SampleRate);
    P = fftshift(P,1);

    TargetFrames = repmat({'0'},FramesPerCycle,1);
    TargetDopplers = zeros(FramesPerCycle,1);
    for Index = 1:size(RangeBins,1)
        TargetFrames{floor(RangeBins(Index,CycleIndex))} = '1';
        TargetFrames{ceil(RangeBins(Index,CycleIndex))} = '1';
        TargetDopplers(floor(RangeBins(Index,CycleIndex))) = ...
            Dopplers(Index,CycleIndex);
        TargetDopplers(ceil(RangeBins(Index,CycleIndex))) = ...
            Dopplers(Index,CycleIndex);
    end
    
    for FrameIndex = 1:FramesPerCycle
        fprintf(FID,'%s,%.1f,%d,%d,%.1f,%.1f', ...
                TargetFrames{FrameIndex}, ...
                TargetDopplers(FrameIndex), ...
                CycleIndex, ...
                FrameIndex, ...
                FrameIndex*MetersPerFrame, ...
                Speed);
        fprintf(FID,',%.4e',10*log10(P(:,FrameIndex)));
        fprintf(FID,'\n');
    end
    
    % Now, let's plot it
    if (Plot)
      useNamedFigure('QuietDataSpectrogram'); clf;
      imagesc(Times,Frequencies,10*log10(P));
      xlabel('Time (sec)');
      ylabel('Frequency (Hz)');
      title(sprintf('Quiet Cycle %d',CycleIndex));
      FileName = sprintf('Figures/QuietCycle-%02d.png',CycleIndex);
    
      print('-dpng',FileName);
    end
end

fclose(FID);
