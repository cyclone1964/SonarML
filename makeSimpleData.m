% So this function makes a simple training data set, using the
% "generateMap" function to make the actual signal. This script
% sets the target strength very high and has no boundary
% reverberation in it to make the investigation easy to visualize
Plot = true;

% The first thing to do is reset the random number generator
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
% target/notarget classes. Let's place targets approximately every
% other frame
TargetsPerCycle = floor(FramesPerCycle/2);

% Now, let's make parameters for each of those cycles. We have the
% sonar speed uniform between -3 and 3, and the Dopplers uniform
% between 10 and 10
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
FID = fopen('SimpleData.csv','w');
fprintf(FID,'Target,Doppler,Cycle,Frame,Range,Speed');
for Index = 1:NumBins
    fprintf(FID,',Bin%02d',Index);
end
fprintf(FID,'\n');

% Now lets go through every cycle and generate the signals, make the
% range doppler map, and write it out.
for CycleIndex = 1:length(Ranges)

    Speed = Speeds(CycleIndex);
    Doppler = Dopplers(CycleIndex);

    % This sets the positions of the targets randomly in bearing at the
    % given ranges.
    Angles = 2*pi*randn(1,TargetsPerCycle);
    Positions = [cos(Angles);sin(Angles);zeros(size(Angles))];
    Positions = repmat(Ranges(:,CycleIndex)',3,1) .* Positions;

    % Generate the signal
    [Samples,Properties] = ...
        generateSamples('Band',Band, ...
                        'PulseLength',PulseLength, ...
                        'TargetStrengths',15*ones(size(Positions)), ...
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
      title(sprintf('Simple Cycle %d',CycleIndex));
      FileName = sprintf('Figures/SimpleCycle-%02d.png',CycleIndex);
    
      print('-dpng',FileName);
    end
end

fclose(FID);
