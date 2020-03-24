%% extract wake and REM sleep

filter_theta = 'n';
review_spikes = 'n';
threshold = 'n';

%% Select Stage Scored File:
working_dir=pwd;
current_dir='C:\SleepData';
cd(current_dir);
scoredCheck = 0;
while isequal(scoredCheck, 0)
    [scoredFile, scoredPath] = uigetfile({'*.xls','Excel 1997-2003 File (*.xls)'},...
        'Select the Sleep Scored File');
    if isequal(scoredFile,0) || isequal(scoredPath,0)
        uiwait(errordlg('You need to select a file. Please try again',...
            'ERROR','modal'));
    else
        cd(working_dir);
        stageScoredFile= fullfile(scoredPath, scoredFile);
        %Load sleep scored file:
        try
            [numData, stringData] = xlsread(stageScoredFile);
            scoredCheck = 1;
        catch %#ok<*CTCH>
            % If file fails to load, it will notify user and prompt to
            % choose another file.
            uiwait(errordlg('Check if the scored file is saved in Microsoft Excel format.',...
             'ERROR','modal'));
         scoredCheck = 0;
        end

    end
end

%% Detect if states are in number or 2-letter format:
if isequal(size(numData,2),3)
    scoredStates = numData(:,2:3);
    clear numData stringData
else
    scoredStates = numData(:,2);
    clear numData
    stringData = stringData(3:end,3);
    [stateNumber] = stateLetter2NumberConverter(stringData);
    scoredStates = [scoredStates stateNumber];
    clear stateNumber stringData
end
epochInSeconds = scoredStates(2,1) - scoredStates(1,1);
startTime = scoredStates(1,1) * 10^6;
endTime = (scoredStates(end,1) + epochInSeconds) * 10^6;

%% Select CSC file:
[CSCFilename, CSCFilePath] = uigetfile({'*.ncs',...
        'Pick CSC files.'},'Select Continuously Sampled Channel File');
cscFile = fullfile(CSCFilePath, CSCFilename);
[TimeStamps, SampleFrequencies, Samples, Header] = Nlx2MatCSC(cscFile, [1 0 1 0 1], 1, 4, [startTime endTime] );

%% Set constants:
Fs = SampleFrequencies(1);
clear SampleFrequencies

% Find AD bit volts in Header:
targ= strfind(Header,'-ADBitVolts');
for i=1:length(targ)
    targIdx(i)= isempty(targ{i}); %#ok<AGROW>
end
ADBitVoltsIdx = find(targIdx==0);   
clear targ targIdx
ADBitVolts = str2double(strrep(Header{ADBitVoltsIdx,1}, '-ADBitVolts ', '')); %#ok<FNDSB>
clear Header

[nsamp,eelen]=size(Samples);

%Reshape the LFP into a single column vector:
newM = nsamp*eelen;
Samples = reshape(Samples, newM, 1);
LFP = Samples * ADBitVolts *10^6; %Convert to microvolts

%normalize signal
% LFP = zscore(LFP);
clear ADBit2uV Samples

%Interpolate time stamps:
interpTimestamps = zeros(eelen*nsamp, 1);
idx = 1;
for i = 1:eelen
  if i < eelen
    t1 = TimeStamps(i);
    t2 = TimeStamps(i+1);
    interval = (t2-t1)/nsamp;
    trange =t1 : interval : t2-interval;
    interpTimestamps(idx:idx+nsamp-1,1) = trange;
  else
    t1 = TimeStamps(i);
    t2 = t1+interval*nsamp;
    trange =(t1 :interval : t2-interval);
    interpTimestamps(idx:idx+nsamp-1,1) = trange;
  end
  idx = idx + nsamp;
end
clear TimeStamps

% Convert from usec to seconds:
timeStamps = interpTimestamps/1000000;
clear interpTimestamps

% Assign states to LFP - output of interest is sleepsamp
lengthSignal = length(LFP);
sleepsamp = zeros(lengthSignal,1);
for i = 1:size(scoredStates, 1)
    if isequal(i, size(scoredStates, 1))
        subIntervalIndx = find(timeStamps >= scoredStates(i,1) & timeStamps < (scoredStates(i,1) + 10));
    else
        subIntervalIndx = find(timeStamps >= scoredStates(i,1) & timeStamps < scoredStates(i+1,1));
    end
    if ~isempty(subIntervalIndx)
        sleepsamp(subIntervalIndx) = scoredStates(i,2);
        clear subIntervalIndx
    end
end


%% LFP and sleepsamp into same variable - note LFP has been zscored

scoredLFP = [LFP,sleepsamp];

i = scoredLFP(:,2);    
data.sleep_idx = i==3;
data.wake_idx = i==1;
data.LFP = scoredLFP(:,1);
data.Fs_LFP = Fs

clear i ADBitVolts ADBitVoltsIdx eelen startTime endTime epochInSeconds t1 t2 trange interval idx
%% filter for theta 

fpass = [1,300];         % input bandpass filter parameters - [1,300] same as Poe 2000
fpass_t = [1,30];
Fs = 1017;
lfp = data.LFP;
N = size(lfp);
dt = 1/Fs;
% time = ((1:N)/Fs)' - dt;


    % filter 1-300Hz
    [b,a] = butter(2,fpass(1)/(Fs/2),'high');
    data.LFP = filtfilt(b,a,lfp);
    % lowpass
    [b,a] = butter(4,fpass(2)/(Fs/2),'low');
    data.LFP = filtfilt(b,a,data.LFP);

if filter_theta == 'y'
    % filter for theta
    [b,a] = butter(2,fpass_t(1)/(Fs/2),'high');
    lfp_theta = filtfilt(b,a,lfp);
    % lowpass
    [b,a] = butter(4,fpass_t(2)/(Fs/2),'low');
    lfp_theta = filtfilt(b,a,lfp_theta);
end

clear LFP N newM nsamp scoredCheck scoredStates b a dt fpass fpass_t 

%% Load cells 

fn = FindFiles('*.ntt','CheckSubdirs', 0);      %prevents code from pulling subfolders in 'cut'
nttNames = cell(size(fn));      %number of ntt files    
ttNum = NaN(size(fn));
for iF = 1:length(fn);
    [pathstr,nttNames{iF}] = fileparts(fn{iF});
    ttNum(iF) = str2num(cell2mat(extractBetween(nttNames{iF},'TT','-')));   %evaluate ntt files to differentiate hippocampal and LC
end
if ~isempty(fn(ttNum<=8))
    S_hc = replay_LoadSpikes('fn',fn(ttNum<=8));  %ntt less than 8 hippocampal - must check with lab notebook
end

clear i iF iC

%% get REM times

detect_theta

run_detect_theta_jk


% Where is zero phase?
[filteredsignal,ts] = filter_data(data.LFP,Fs,'bandpass',3,5,10);
phase0 = abs(IP_theta.data)<.01;
sum(phase0);

 figure; plot(timeStamps,filteredsignal,'linewidth',1); hold on; plot(timeStamps(phase0),filteredsignal(phase0),'r.','markersize',15)



% restrict spkTimes to REM periods    
    S_hc_REM = S_hc;    %spikes have been loaded for hippocampus
    for iC = 1:length(S_hc);    %restrict spikes to only run times (excludes between runs and sleep)
        S_hc_REM{iC} = S_hc_REM{iC}.restrict(REM_on,REM_off);   %new spikes specific to runs variable created
    end

    % restrict spkTimes to wake periods
    S_hc_wake = S_hc;    %spikes have been loaded for hippocampus
    for iC = 1:length(S_hc);    %restrict spikes to only run times (excludes between runs and sleep)
        S_hc_wake{iC} = S_hc_wake{iC}.restrict(wake_on,wake_off);   %new spikes specific to runs variable created
    end

%% plot detected spikes 

if review_spikes == 'y'
    spktimes_REM = S_hc_REM{1}.data;
    theta_plot = REM_theta_trough.data;
    wakeREM.theta_troughs = ismember(timeStamps,theta_plot);

    cc = get(gca,'ColorOrder');

for i=1:length(spktimes_REM);
           
        clf
    
        tmin = spktimes_REM(i)-5;
        tmax = spktimes_REM(i)+5;
        
        axh = gca;
        plot(timeStamps,wakeREM.LFP)
        hold on
        line([spktimes_REM(i) spktimes_REM(i)], get(axh,'YLim'), 'Color', 'm', 'LineStyle','-');
        
        
        hold on
        plot(timeStamps(wakeREM.REM),repmat(1*max(wakeREM.LFP),sum(wakeREM.REM),1),'.','Color',cc(3,:))
        hold on
        plot(timeStamps(wakeREM.theta_troughs),repmat(1*min(wakeREM.LFP),sum(wakeREM.theta_troughs),1),'k*')
        axis([tmin tmax  get(axh,'YLim') get(axh,'YLim')]);
    
    
    pause
end
end

