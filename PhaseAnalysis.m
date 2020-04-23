%% Created by R.G. 2020 for the purpose of analyzing spike-phase preference 
%% during active behavior and subsequent sleep 

%% Description of code:
    %   1. Load the Poe Lab sleep scored file for the animal/day of
    %   interest
    %
    %   2. Load the CSC file to be used in phase analysis 
    %
    %   3. Package variables into 'data' structure containing - 
    %        sleep_idx = logical containing REM epochs [REM=1;wake=0] 
    %        wake_idx = local containing wake epochs [wake=1]
    %        LFP = LFP in uV
    %        Fs_LFP = Sampling rate of recording
    %
    %   4. Bandpass the signal from fpass = [1,300]Hz
    %
    %   5. Get the Hilbert transform and package into tsd format with -  
    %       IP = instantaneous phase (radians)
    %       IA = instantaneous amplitude (power)
    %       IF = instantaneous frequency 
    %      Zscore IA
    %
    %   6. Load spikes into tsd format from .ntt file
    %
    %   7. Obtain maze run and REM timestamps 
    %
    %   8. Restrict spikes to maze run times 
    %
    %   9. Restrict spikes to REM periods
    %
    %   10. Obtain a threshold for theta amplitude using only top 85th
    %   percentile 
    %
    %   11. Restrict spikes by theta amplitude and convert from tsd to cell
    %
    %   12. Package hilbert and set parameters for phase analysis
    %
    %   13. Using circular statistics toolbox - 
    %       Bin spikes during waking run of interest (familiar or altered)
    %       by theta phase and obtain polar plot of spike phase.
    %       Bin spikes during REM sleep by theta phase and obtain polar
    %       plot of spike phase. 
    %
    %   14. Circstat outputs of interest -
    %       mu    Mean phase angle
    %       pval  rayleigh test for significance of phase preference 
    %       r     vector length    
                

%% Determine auxilary outputs
filter_theta = 'n';
review_spikes = 'n';
threshold = 'n';
plotphase = 'n';

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

% Detect if states are in number or 2-letter format:
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
[TimeStamps, SampleFrequencies, Samples, Header] = Nlx2MatCSC(cscFile, [1 0 1 0 1], 1, 1, [startTime endTime] );

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
data.Fs_LFP = Fs;

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



%% Get Hilbert and package amplitude (aribitrary units) and phase (r)

CSC = tsd(timeStamps,data.LFP);
[IF_theta, IA_theta, IP_theta, CSC0_theta]=InstSig_theta(CSC);

theta_amp_zsc = zscore(IA_theta.data);              % get the zscored theta amp
theta_amp = tsd(IA_theta.range,theta_amp_zsc);      % put zscored theta amp back into tsd format

clear theta_amp_zsc
%% Load spikes

fn = FindFiles('*.ntt','CheckSubdirs', 0);      %prevents code from pulling subfolders in 'cut'
nttNames = cell(size(fn));      %number of ntt files    
ttNum = NaN(size(fn));
for iF = 1:length(fn);
    [pathstr,nttNames{iF}] = fileparts(fn{iF});
    ttNum(iF) = str2double(cell2mat(extractBetween(nttNames{iF},'TT','-')));   %evaluate ntt files to differentiate hippocampal and LC
end
if ~isempty(fn(ttNum<=8))
    S_hc = replay_LoadSpikes('fn',fn(ttNum<=8));  %ntt less than 8 hippocampal - must check with lab notebook
end
if ~isempty(fn(ttNum>8));
    S_LC = replay_LoadSpikes('fn',fn(ttNum>8)); %ntt files greater than 8 LC - must check with lab notebook
end

clear ttNum nttNames
%% Generate run and REM times

t = load('SessionTimes.mat'); %loads session times .mat file that was created when the replay_MarkBehaviorTimes.m function creates (saved into host folder)
run_on = [t.session1_times([2,4,6]); t.session2_times([2,4,6])]; %uses only the on track times created in MarkBehavior
run_off = [t.session1_times([3,5,7]); t.session2_times([3,5,7])]; %uses only the off track times created in MarkBehavior

% get familiar and altered maze run separetely

run_on_fam = t.session1_times([2,4,6]);
run_off_fam = t.session1_times([3,5,7]);

run_on_alt = t.session2_times([2,4,6]);
run_off_alt = t.session2_times([3,5,7]);

clear t

% bin REM epochs
% first set parameters
maxTimeBetween = .005;
maxdist = maxTimeBetween * Fs;

timestamps_REM = timeStamps(data.sleep_idx);

REM_on= timestamps_REM(1);
REM_off = [];
for i = 1:length(timestamps_REM) -1;    
     if timestamps_REM(i+1)-timestamps_REM(i) > maxdist;
         REM_on = [REM_on; timestamps_REM(i+1)];
         REM_off = [REM_off; timestamps_REM(i)];
     end
end

REM_on(end) = [];

%% restrict spikes to familiar and altered maze runs
S_hc_fam = S_hc;
S_hc_alt = S_hc;
for iC = 1:length(S_hc);
    S_hc_fam{iC} = S_hc_fam{iC}.restrict(run_on_fam,run_off_fam);
    S_hc_alt{iC} = S_hc_alt{iC}.restrict(run_on_alt,run_off_alt);
end

%% restrict spikes to REM periods
S_hc_REM = S_hc;    %spikes have been loaded for hippocampus
for iC = 1:length(S_hc);    %restrict spikes to only run times (excludes between runs and sleep)
    S_hc_REM{iC} = S_hc_REM{iC}.restrict(REM_on,REM_off)';   %new spikes specific to runs variable created
end

clear S_hc

%% set threshold - work in progress
mnl_parm=[85];
thresh = prctile(theta_amp.data,mnl_parm); 

%% convert tsd to cell - definitely some redundancies 

spikeAmp_wake = [];
spikeTimes_wake = cell(length(S_hc_alt),1);
for i = 1:length(S_hc_alt);
    spikeTimes_wake{i} = [S_hc_alt{i}.data];
    spikeAmp_wake{i} = theta_amp.data(S_hc_alt{i}.data);   % get amplitude at spike times 
    k{i} = spikeAmp_wake{i} > thresh;                           % threshold step
    allspikes = S_hc_alt{i}.data;                        
    S_hc_alt_thresh{i} = ts(allspikes(k{i}));              % put thresholded spiketimes back into tsd format
end

spikeTimes_wake = cell(length(S_hc_alt_thresh),1);
for i = 1:length(S_hc_alt_thresh);
    spikeTimes_wake{i} = [S_hc_alt_thresh{i}.data];
end

clear allspikes k 

spikeAmp_rem = [];
spikeTimes_rem = cell(length(S_hc_REM),1);
for i = 1:length(S_hc_REM);
    spikeTimes_rem{i} = [S_hc_REM{i}.data];
    spikeAmp_rem{i} = theta_amp.data(S_hc_REM{i}.data);   % get amplitude at spike times 
    k{i} = spikeAmp_rem{i} > thresh;                           % threshold step
    allspikes = S_hc_REM{i}.data;                        
    S_hc_rem_thresh{i} = ts(allspikes(k{i}));              % put thresholded spiketimes back into tsd format
end

spikeTimes_rem = cell(length(S_hc_rem_thresh),1);
for i = 1:length(S_hc_rem_thresh);
    spikeTimes_rem{i} = [S_hc_rem_thresh{i}.data]';
end

clear allspikes k S_hc_fam_thresh S_hc_rem_thresh

%% Get phase from tsd hilbert and set parameters for phase analysis
signalPhase = IP_theta.data;
timeStamps_sp = IP_theta.range;

allPhases = -pi:pi/12:pi;
PhaseCenters = allPhases(1:end-1)+pi/24;

phasePi = abs(signalPhase-pi)<0.01;
phase0 = abs(signalPhase)<.01;

%% plot the phase
if plotphase == 'y'

    f = figure;
    axh = gca;
    cc = get(gca,'ColorOrder');

    plot(timeStamps,lfp_theta,'linewidth',1); hold on; plot(timeStamps(phase0),lfp_theta(phase0),'r.','markersize',15) % plot LFP with phase indicated
    hold on
    plot(timeStamps(data.sleep_idx),repmat(1*max(lfp_theta),sum(data.sleep_idx),1),'.','Color',cc(2,:))                % plot REM epochs above LFP
    hold on
    plot(timeStamps(data.wake_idx),repmat(1*max(lfp_theta),sum(data.wake_idx),1),'.','Color',cc(3,:))                  % plot wake epochs above LFP
end

%% Wake circstat results
for i = 1:length(spikeTimes_wake)
    spikePhases_wake{i} = IP_theta.data(spikeTimes_wake{i});
    spikecounts_wake{i}= histcounts(spikePhases_wake{i},allPhases)';
end

for i=1:length(spikeTimes_wake)
    clf
    subplot(1,2,1);bar(PhaseCenters,spikecounts_wake{i});title(sprintf('Cell%d',i));
    subAx(i) = subplot(1,2,2);polar(PhaseCenters,spikecounts_wake{i}'); hold on

    mu_wake(i) = circ_mean(PhaseCenters',spikecounts_wake{i});
    
    %rayleigh test
    p_wake(i) = circ_rtest(PhaseCenters,spikecounts_wake{i});    
    h_wake = polar(subAx(i),[0 mu_wake(i)],[0 max(spikecounts_wake{i})/2]); set(h_wake,'linewidth',2);title(p_wake(i))

    %get vector length
    r_wake_raw(i) = circ_r(PhaseCenters',spikecounts_wake{i});
    
    pause
end

%% REM circstat results 

for i = 1:length(spikeTimes_rem)
    spikePhases_rem{i} = IP_theta.data(spikeTimes_rem{i});
    spikecounts_rem{i}= histcounts(spikePhases_rem{i},allPhases)';
end

for i=1:length(spikeTimes_rem)
    clf
    subplot(1,2,1);bar(PhaseCenters,spikecounts_rem{i});title(sprintf('Cell%d',i));
    subAx(i) = subplot(1,2,2);polar(PhaseCenters,spikecounts_rem{i}'); hold on

    mu_rem(i) = circ_mean(PhaseCenters',spikecounts_rem{i});
    
    %rayleigh test
    p_rem(i) = circ_rtest(PhaseCenters,spikecounts_rem{i});    
    h_rem = polar(subAx(i),[0 mu_rem(i)],[0 max(spikecounts_rem{i})/2]); set(h_rem,'linewidth',2);title(p_rem(i))
    
    %get vector length
    r_rem_raw(i) = circ_r(PhaseCenters',spikecounts_rem{i});
    
    pause
end

%% workshop - where I put code to be later incorporated into script

