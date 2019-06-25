% This script makes an advanced data set. This data set has 
% 128 pt spectra, a pulse length commensurate with that length,
% targets NOT aligned with frames, boundary reverberation enabled
% with a beam pattern corresponding to an aperture of 2
% wavelengths and target ranges. 
%
% THe parameters of this data set are chosen to match a known
% in-water data set (at least to start), which has a target limited
% to 4 knots and a pulse length of 64 ms. So we set the band to be
% 1k wide and the frame size to be 64. 
Plot = 0;

% The first thing to do is reset the random number generator
rng('default');

% This is the number of frequency bins in a frame
NumBins = 64;

% For this simulation, the band is between 18 and 22 k
Band = [19500 20500];

% The sound speed is about this
SoundSpeed = 1473;

% That means that the complex baseband sample rate is the bandwidth
SampleRate = diff(Band);
CenterFrequency = mean(Band);
PulseLength = NumBins/SampleRate;
SecondsPerFrame = NumBins/SampleRate;
MetersPerFrame = SecondsPerFrame * SoundSpeed/2;

% We want cycles that are 2 seconds long, so it's this many samples,
% but we round up to an integer number of frames
SamplesPerCycle = 2*SampleRate;
FramesPerCycle = ceil(SamplesPerCycle/NumBins);
SamplesPerCycle = FramesPerCycle*NumBins;
SecondsPerCycle = FramesPerCycle * SecondsPerFrame;
MetersPerCycle = FramesPerCycle*MetersPerFrame;
FrameRanges = (1:FramesPerCycle)*MetersPerFrame;

HzPerBin = SampleRate/NumBins;
MpsPerBin = HzPerBin * SoundSpeed/CenterFrequency;

% now make some pictures for illustrative purposes
if (Plot > 1)
    WaterColumn = initializeWaterColumn('SoundSpeed',SoundSpeed, ...
                                        'ScatteringStrength',-75);
    SurfaceBoundary = ...
        initializeSurfaceBoundary('SeaState',3,'WaveHeight',3,'WindSpeed',10);
    BottomBoundary = ...
        initializeBottomBoundary('Depth',600,'Porosity',0.5,'GrainSize',5);
    Environment = initializeEnvironment('SurfaceBoundary',SurfaceBoundary, ...
                                        'BottomBoundary',BottomBoundary, ...
                                        'WaterColumn',WaterColumn);
    PlatformDepth = 200;
    
    BeamParameters = ...
        {[], [], 'No Beam'
         4, 10, 'Full Aperture Beam'};
    ComponentParameters = {'Volume' [0] [-100; -100]
                        'Surface' [-100] [0; -100]
                        'Bottom' [-100] [-100; 0]
                        'Total' [0] [0; 0]};
    
    for BeamIndex = 1:size(BeamParameters,1)
        useNamedFigure('Level'); clf; hold on;
        for ComponentIndex = 1:size(ComponentParameters,1)
            rng('default');
            [Samples,Properties] = ...
                generateSamples('Band',Band, ...
                                'PulseLength',PulseLength, ...
                                'Aperture',BeamParameters{BeamIndex,1}, ...
                                'Baffling',BeamParameters{BeamIndex,2}, ...
                                'Steering',[1;0;0], ...
                                'Environment',Environment, ...
                                'CycleLength',SecondsPerCycle, ...
                                'PlatformSpeed',10, ...
                                'PlatformDepth',PlatformDepth, ...
                                'VolumeReverbAdjustment', ...
                                ComponentParameters{ComponentIndex,2}, ...
                                'BoundaryReverbAdjustment', ...
                                ComponentParameters{ComponentIndex,3});
            Times = (1:length(Samples))/diff(Band);
            Temp = conv(abs(Samples),ones(100,1)/100,'same');
            plot(Times,20*log10(Temp));
        end
        Temp = axis; Temp(3:4) = [100 180]; axis(Temp);
        xlabel('Time (s)'); ylabel('Level (dB)');
        legend('Volume','Surface','Bottom','Total');
        
        BeamName = BeamParameters{BeamIndex,3};
        title(['Improved Model: ' BeamName])
        prettyPlot;
        print('-dpng', ...
              ['ImprovedModelEnvelope-',strrep(BeamName,' ','') '.png']);
    
        useNamedFigure('NomimalBeam'); clf; 
        [S,Frequencies,Times,P] = ...
            spectrogram(Samples,NumBins,0,NumBins,SampleRate);
        P = fftshift(P,1);
        Doppler = SoundSpeed*(Frequencies-mean(Frequencies))/CenterFrequency;
        imagesc(Times,Doppler,10*log10(P));
        set(gca,'YDir','normal');
        title(['Improved Model: ' BeamName]);
        print('-dpng',...
              ['ImprovedModelSpectrogram-' strrep(BeamName,' ','') '.png']);
    end
end

% Now an environment meant to mimic the target situation
WaterColumn = initializeWaterColumn('SoundSpeed',SoundSpeed, ...
                                    'ScatteringStrength',-75);
SurfaceBoundary = ...
    initializeSurfaceBoundary('SeaState',3,'WaveHeight',3,'WindSpeed',5);
BottomBoundary = ...
    initializeBottomBoundary('Depth',400,'Porosity',0.5,'GrainSize',4);
DeepEnvironment = initializeEnvironment('Surface',SurfaceBoundary, ...
                                        'Bottom',BottomBoundary, ...
                                        'WaterColumn',WaterColumn);
ShallowEnvironment = DeepEnvironment;
ShallowEnvironment.Bottom.Depth = 150;

% Now, the number of Cyles. We want a total of 16k training
% spectra. 
NumCycles = ceil((8*1024)/FramesPerCycle);

DataSetParameters = ...
    {'Quiet', [], [], [-100; -100] 
     'Unaligned', [], [], [-100; -100]
     'Boundary', [], [], [0;0]
     'BeamPattern', 4, 10, [0; 0],
     'BeamShallow', 4, 10, [0; 0]};

for DataSetIndex = 1:size(DataSetParameters,1)
    
    % Extract the data set information
    DataSetName = DataSetParameters{DataSetIndex,1};
    Aperture = DataSetParameters{DataSetIndex,2};
    Baffling = DataSetParameters{DataSetIndex,3};
    BoundaryAdjustement = DataSetParameters{DataSetIndex,4};
    fprintf('Process %s \n',DataSetName);

    if (DataSetIndex == size(DataSetParameters,1))
        Environment = ShallowEnvironment;
    else
        Environment= DeepEnvironment;
    end
    
    % Open the file and write the header
    rng('default');
    FID = fopen([DataSetName 'Data.csv'],'w');

    % Now lets go through every cycle and generate the signals, make the
    % range doppler map, and write it out. Let's also save the images
    % of the cycles.
    for CycleIndex = 1:NumCycles

        % Get a set of ranges for each target. These are uniformly
        % centered around a set of range bin centers to keep them from
        % being at the same range.
        FrameIndices = 3:4:(FramesPerCycle-2);
        if (DataSetIndex == 1) 
            FrameIndices = FrameIndices + ...
                randi(3,size(FrameIndices)) - 2;
        else
            FrameIndices = (FrameIndices + (2*rand(size(FrameIndices)) - 1));
        end
        Ranges = MetersPerFrame * FrameIndices;

        % Set the platform speed and the dopplers 
        PlatformSpeed = 5*rand(1);
        Bearings = 360*rand(size(Ranges));
        Velocities = ...
          20*[cosd(Bearings)
              sind(Bearings)
              zeros(size(Bearings))];

        % Now, get a random direction that the beam will be steered in.
        Azimuth = (60*rand(1)-30) * pi/180;
        Elevation = (20*rand(1)-10) * pi/180;
        Direction = computeDirection([0 Elevation Azimuth]');
        Positions = Direction * Ranges;

        % Set the water depth and platform depth
        PlatformDepth = 50+rand*(Environment.Bottom.Depth-100);
        [Samples,Properties] = ...
            generateSamples('Band',Band, ...
                            'PulseLength',PulseLength, ...
                            'Aperture',Aperture, ...
                            'Baffling',Baffling, ...
                            'Steering',Direction, ...
                            'Environment',Environment, ...
                            'CycleLength',SecondsPerCycle, ...
                            'TargetStrengths',-5*ones(size(Ranges)), ...
                            'TargetPositions',Positions, ...
                            'TargetVelocities',Velocities, ...
                            'PlatformSpeed',PlatformSpeed, ...
                            'PlatformDepth',PlatformDepth, ...
                            'VolumeReverbAdjustment',0, ...
                            'BoundaryReverbAdjustment',BoundaryAdjustement);
    
        % Make a spectrogram
        [S,Frequencies,Times,P] = ...
            spectrogram(Samples,NumBins,0,NumBins,SampleRate);
        P = fftshift(P,1);

        TargetFrames = repmat({'0'},FramesPerCycle,1);
        TargetDopplers = Direction' * Velocities;
        SourceDoppler = PlatformSpeed * Direction(1);
        TotalDopplers = TargetDopplers + SourceDoppler;
        TargetBins = round(TotalDopplers/MpsPerBin) + NumBins/2 + 1;
        TargetDopplers = zeros(size(Times));
        for Index = 1:length(Ranges)
          if (DataSetIndex == 1)
            FrameIndex = round(Ranges(Index)/MetersPerFrame)+1;
            TargetFrames{FrameIndex} = '1';
            TargetDopplers(FrameIndex) = TargetBins(Index);
          else
            FrameIndex = floor(Ranges(Index)/MetersPerFrame)+1;
            TargetFrames{FrameIndex} = '1';
            TargetDopplers(FrameIndex) = TargetBins(Index);
            FrameIndex = ceil(Ranges(Index)/MetersPerFrame)+1;
            
            TargetFrames{FrameIndex} = '1';
            TargetDopplers(FrameIndex) = TargetBins(Index);
          end
        end
        
        for FrameIndex = 1:FramesPerCycle
            fprintf(FID,'%s,%.1f,%d,%d,%.1f,%.1f,%.0f,%.0f', ...
                    TargetFrames{FrameIndex}, ...
                    TargetDopplers(FrameIndex), ...
                    CycleIndex, ...
                    FrameIndex, ...
                    FrameRanges(FrameIndex), ...
                    PlatformSpeed, ...
                    Azimuth * 180/pi, ...
                    Elevation * 180/pi);
            fprintf(FID,',%.4e',10*log10(P(:,FrameIndex)));
            fprintf(FID,'\n');
        end
        
        % Now, let's plot it
        if ((Plot > 0  && CycleIndex < 10) || CycleIndex == 1)
            useNamedFigure([DataSetName 'Spectrogram']); clf;
            Doppler = ...
                SoundSpeed*(Frequencies-mean(Frequencies))/CenterFrequency;
            imagesc(Times,Doppler,10*log10(P));
            xlabel('Time (sec)');
            ylabel('Frequency (Hz)');
            set(gca,'YDir','normal');
            title(sprintf('%s Cycle %d - %.1f mps, z = %.0f, d = %.0f', ...
                          DataSetName, ...
                          CycleIndex, ...
                          PlatformSpeed, PlatformDepth, ...
                          Environment.Bottom.Depth));
            FileName = sprintf('Figures/%sCycle-%02d.png', ...
                               DataSetName,CycleIndex);
            prettyPlot(); print('-dpng',FileName);

            if (Plot > 1 || CycleIndex == 1)
                useNamedFigure([DataSetName 'Envelope']); clf
                plot(Times,10*log10(sum(P,1)));
                xlabel('Time (sec)'); ylabel('Level (dB)');
                set(gca,'YLim',[60 150]);
                title(sprintf('%s Cycle %d - %.1f mps, z = %.0f, d = %.0f', ...
                              DataSetName, ...
                              CycleIndex, ...
                              PlatformSpeed, PlatformDepth, ...
                              Environment.Bottom.Depth));
                FileName = sprintf('Figures/%sEnvelope-%02d.png', ...
                                   DataSetName, CycleIndex);
                prettyPlot(); print('-dpng',FileName);          
            end
        end
    end
end

fclose(FID);
