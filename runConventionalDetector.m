% This script runs a "conventional" detector on all the cycles in the
% data sets. This involves computing a normalizer on those cycles
% using an inverse template, normalzing the cycle, and then finding
% threshold crossings. 

% The sample rate as defined in the simulation and the list of 
% data set names
SampleRate = 1000;
DataSetNames = {'QuietData', ...
             'UnalignedData', ...
             'BoundaryData', ...
             'BeamPatternData', ...
             'BeamShallowData'};

% This is the threshold, in dB that seems to reach the best balance
% of FA/MD
Threshold = 6;

% Do each data set name
for DataSetName = DataSetNames
  
  % Load the data: all of it!
  DataSetName = DataSetName{1};
  Data = load([DataSetName '.csv']);

  % Compute the number of cycles and the number of frames per
  % cycle, and some other stuff
  NumCycles = length(unique(Data(:,3)));
  FramesPerCycle = length(unique(Data(:,4)));
  BinsPerFrame = size(Data,2) - 8;
  SecPerFrame = BinsPerFrame/SampleRate;
  HzPerBin = SampleRate/BinsPerFrame;

  % Define the template: I tried lots of parameters for this and
  % this seems to be the best.
  Template = ones(5,9);
  Template(2:4,3:7) = 0;
  Template = Template /sum(Template(:));
  StartingIndex = 1;

  % Now compute the number of correct classifications, false
  % alarms, and missed detections by processing each cycle
  NumCorrect = 0;NumFA = 0;NumMD = 0;
  for CycleIndex = 1:NumCycles

    % This is the indices for this particular cycle
    Indices = (1:FramesPerCycle) + StartingIndex - 1;
    StartingIndex = StartingIndex + FramesPerCycle;

    % Get the cycle bins
    CycleBins = 10.^(Data(Indices,9:end)/10);
    
    % Augment front, back, left, and right to pad it out to get the
    % normalization proper
    PaddedBins = [repmat(CycleBins(1,:),2,1)
                  CycleBins
                  repmat(CycleBins(end,:),2,1)];
    PaddedBins = [repmat(PaddedBins(:,1),1,4) ...
                  PaddedBins ...
                  repmat(PaddedBins(:,end),1,4)];

    % Compute the normalizer by convolving the padded bins with the
    % template, then re-sizing to select the original ones
    Normalizer = conv2(PaddedBins,Template,'same');
    Normalizer = Normalizer(3:end-2,5:end-4);
    
    % Comput the SNR and convert to dB
    SNR = CycleBins ./ Normalizer;
    SNR = 10*log10(SNR);

    % WE use the time and frequencies in the plotting blow
    Times = (1:FramesPerCycle) * SecPerFrame;
    Frequencies = (1:BinsPerFrame) * HzPerBin;
    Frequencies = Frequencies - mean(Frequencies);

    % Let's plot just the first one and print it to a file
    if (CycleIndex == 1)
      useNamedFigure('DataSetName'); clf;
      imagesc(Times, Frequencies, SNR');
      set(gca,'YDir','reverse');
      xlabel('Time (sec)');
      ylabel('Frequency (Hz)'); 
      title(sprintf('SNR: %s',DataSetName));
      prettyPlot
      print('-dpng',['Figures/' DataSetName '-SNR.png']);
    end

    % Find those entries that have a SNR over the threshold, recuce
    % to only those frames with detections in them
    [Frames, Bins] = find(SNR > Threshold);
    Frames = unique(Frames);

    % This is a vector that has ones in it where there are
    % detections, and is thus comparable to the labels in the first
    % column of the data.
    Prediction = zeros(FramesPerCycle,1);
    Prediction(Frames) = 1;

    % Update the statistics
    NumCorrect = NumCorrect + ...
        length(find(Prediction == Data(Indices,1)));
    NumFA = NumFA + length(find(Prediction ==1 & Data(Indices,1) == 0));
    NumMD = NumMD + length(find(Prediction ==0 & Data(Indices,1) ==1));
    
  end

  % Print the results
  fprintf('File %s Thres %.0f = %.2f (FA: %.0f,MD: %.0f)\n', ...
          DataSetName,  ...
          Threshold, ...
          100*NumCorrect/size(Data,1), ...
          100*NumFA/size(Data,1),  ...
          100*NumMD/size(Data,1));
end