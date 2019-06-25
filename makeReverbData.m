% So this function makes a simple training data set. It makes
% active sonar cycles using a "filtered noise" model. This model
% makes STFT frames of a specific length, multiplies them by
% complex random noise, IFFTs it, then blends the frames together
%
% In this model, the transmit beam pattern is fully omni, the
% platform is moving, as is the target. 
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

% We want cycles that are 1 second long, so it's this many samples
SamplesPerCycle = SampleRate;
FramesPerCycle = ceil(SamplesPerCycle/NumBins);
SamplesPerCycle = FramesPerCycle*NumBins;

% Now, the number of Cyles. We want a total of 16k training
% spectra. 
NumCycles = ceil((8*1024)/FramesPerCycle);

% In order for the trianing to work, we want multiple targets per
% cycle. Otherwise there is a 32 to 1 imbalance in the
% target/notarget classes. So we place 8 targets per cycle
TargetsPerCycle = 64;

% Now, let's make parameters for each of those cycles. We have the
% sonar speed uniform between -5 and 5, the target range uniform
% between 100 and 600, and the Dopplers uniform between -21 and 21
Speeds = 6*randn(TargetsPerCycle,NumCycles)-3;
Dopplers = 20 * rand(TargetsPerCycle,NumCycles)-10;
WaterDepth = 250;
PlatformDepths = WaterDepth * rand(NumCycles);
MetersPerFrame = 750*SecondsPerFrame;

% Now, we need the range fixed to specific frames, so let's do that
% now
RangeBins = zeros(size(Dopplers));
for CycleIndex =1:NumCycles
    RangeBins(:,CycleIndex) = randperm(FramesPerCycle-4, TargetsPerCycle)+2;
end
Ranges = MetersPerFrame * (RangeBins-1);

% Open the file and write the header
FID = fopen('AdvancedData.csv','w');
fprintf(FID,'Target,Doppler,Cycle,Frame,Range,Speed');
for Index = 1:NumBins
    fprintf(FID,',Bin%02d',Index);
end
fprintf(FID,'\n');

% Now lets go throug every cycle and genrate the signals, make the
% range doppler map, and write it out. Let's also save the images
% of the cycles.
for CycleIndex = 1:length(Ranges)

    Speed = Speeds(CycleIndex);
    Doppler = Dopplers(CycleIndex);
    
    Angles = 2*pi*randn(1,TargetsPerCycle);
    Positions = [cos(Angles);sin(Angles);zeros(size(Angles))];
    Positions = repmat(Ranges(:,CycleIndex)',3,1) .* Positions;
    
    [Samples,Properties] = ...
        generateMap('Band',Band, ...
                    'WaterDepth',WaterDepth, ...
                    'PlatformDepth',PlatformDepths(CycleIndex), ...
                    'PulseLength',PulseLength, ...
                    'TargetStrengths',20*ones(size(Positions)), ...
                    'TargetPositions',Positions, ...
                    'TargetDopplers',Dopplers(:,CycleIndex), ...
                    'PlatformSpeed',Speed);
    
    % Make a spectrogram
    SoundSpeed = Properties.SoundSpeed;

    [S,Frequencies,Times,P] = spectrogram(Samples,32,0,32,SampleRate);
    P = fftshift(P,1);

    TargetFrames = repmat({'NoTarget'},FramesPerCycle,1);
    TargetDopplers = zeros(FramesPerCycle,1);
    for Index = 1:size(RangeBins,1)
        TargetFrames{RangeBins(Index,CycleIndex)} = 'Target';
        TargetDopplers(RangeBins(Index,CycleIndex)) = ...
            Dopplers(Index,CycleIndex);
    end
    
    for FrameIndex = 1:FramesPerCycle
        fprintf(FID,'%s,%.1f,%d,%d,%.1f,%.1f', ...
                TargetFrames{FrameIndex}, ...
                TargetDopplers(FrameIndex), ...
                CycleIndex, ...
                FrameIndex, ...
                Ranges(FrameIndex), ...
                Speed);
        fprintf(FID,',%.4e',P(:,FrameIndex));
        fprintf(FID,'\n');
    end
    
    % Now, let's plot it
    if (Plot)
      useNamedFigure('SimpleDataSpectrogram'); clf;
      imagesc(Times,Frequencies,10*log10(P));
      xlabel('Time (sec)');
      ylabel('Frequency (Hz)');
      title(sprintf('Training Cycle %d',CycleIndex));
      FileName = sprintf('TrainingCycle-%02d.png',CycleIndex);
    
      print('-dpng',FileName);
    end
end

fclose(FID);
