% close all
BCI2kPath = '/Users/nkb/Documents/NCAN/BCI2000tools';
addpath('/Users/nkb/Documents/NCAN/code/MATLAB_tools')
bci2ktools(BCI2kPath);


sessionID = "TI_PNS_20240612";
sessionID = "TI_PNS_20240625";

dataPath = fullfile('/Users/nkb/Library/CloudStorage/Box-Box/Temporal-Interference-Project/data',sessionID);
metapath = fullfile("/Users/nkb/Library/CloudStorage/Box-Box/Temporal-Interference-Project/data",sessionID,"PNS_TI_Runsheet.xlsx");
% metapath = "/Users/nkb/Library/CloudStorage/Box-Box/Temporal-Interference-Project/data/TI_PNS_20240507/recording_20240507_NBMP.csv";
runMetaData = readtable(metapath,'VariableNamingRule','preserve');
runMetaData = table2struct(runMetaData);
datFiles = dir(fullfile(dataPath,'*.dat'));
dataStruct = struct();
meta_out = struct();
EMG_IDS = {'flex','ext','thumb'}; % add to this if more EMG channels are parsed. 
for i=1:length(datFiles)
    processed = struct();
    [sig,state,param] = load_bcidat(fullfile(datFiles(i).folder,runMetaData(i).fname),'-calibrated');
    dataStruct(i).raw = sig;
    dataStruct(i).state = state;
    dataStruct(i).param = param;
    dataStruct(i).fs = param.SamplingRate.NumericValue;
    dataStruct(i).sampleBlock = param.SampleBlockSize.NumericValue;
    dataStruct(i).stimType = runMetaData(i).StimType;
    dataStruct(i).deltaF = runMetaData(i).deltaF;
    chans = param.ChannelNames.Value;
    containsSubstring = any(contains(chans, EMG_IDS), 2);
    emg_chan = chans(containsSubstring);
    eeg_chan = chans(~containsSubstring);
    emgSigs = sig(:,containsSubstring);
    EEGSigs = sig(:,~containsSubstring);
    C = unique(regexprep(emg_chan,'[^a-zA-Z]',''),'stable');
    [raw,filt] = process_EMG(emgSigs,C,dataStruct(i).fs);
    dataStruct(i).signals = raw;
    dataStruct(i).signals_f = filt;
    % TODO: process EEG 
    dataStruct(i).EEG = reref_EEG(EEGSigs,eeg_chan);
    

end
%%
figure(1)
rows = length(dataStruct);
cols = length(dataStruct(1).signals);
count = 1;
for i=1:rows
    for j = 1:cols
        subplot(rows,cols, count);set(gca,'ButtonDownFcn',@fig_from_subplot);
        count =count+1;
        t = linspace(0,length(dataStruct(i).signals(j).data)/dataStruct(i).fs,length(dataStruct(i).signals(j).data));
        % plot(t,dataStruct(i).signals(j).data)
        hold on
        plot(t,dataStruct(i).signals_f(j).data,'red');
        hold off
        stdev = std(dataStruct(i).signals_f(j).data);
        ylim([-4*stdev, 4*stdev]);
        xlabel('time (s)');
        ylabel('amplitude (uV)')
        metadata = sprintf("%s Run %d, %s Stim, \n deltaF=%g Hz, amp=%s \n%s",dataStruct(i).signals(j).name,runMetaData(i).Run,runMetaData(i).StimType, runMetaData(i).deltaF, runMetaData(i).Amplitude, runMetaData(i).Notes);
        title(metadata);
    end
end
figure(2)
count = 1;
for i=1:rows
    for j = 1:cols
        subplot(rows,cols, count);set(gca,'ButtonDownFcn',@fig_from_subplot);
        count =count+1;
        t = linspace(0,length(dataStruct(i).signals(j).data)/dataStruct(i).fs,length(dataStruct(i).signals(j).data));
        % plot(t,dataStruct(i).signals(j).data)
        hold on
        plot(t,dataStruct(i).signals(j).data,'blue');
        hold off
        % stdev = std(dataStruct(i).signals_f(j).data);
        % ylim([-4*stdev, 4*stdev]);
        xlabel('time (s)');
        ylabel('amplitude (uV)')
        metadata = sprintf("%s Run %d, %s Stim, \n deltaF=%g Hz, amp=%s \n%s",dataStruct(i).signals(j).name,runMetaData(i).Run,runMetaData(i).StimType, runMetaData(i).deltaF, runMetaData(i).Amplitude, runMetaData(i).Notes);
        title(metadata);
    end
end
%% Peak Detection
for i=1:length(dataStruct)
    sample_offset = 1000;
    data = dataStruct(i).signals(1).data(sample_offset:end-sample_offset);
    X = get_artifact_peaks(abs(data),max(data)*0.75,dataStruct(i).fs,0);
    % dataStruct(i).emg_epochs = zeros(length(dataStruct(i).signals_f));
    for j=1:length(dataStruct(i).signals_f)
        samplen = 2150;
        onset_samps = 500;
        epoch_time = linspace(-(onset_samps/dataStruct(i).fs),(samplen-onset_samps)/dataStruct(i).fs,samplen);
        epochs = epoch_by_peak_loc(dataStruct(i).signals_f(j).data(1000:end-1000),X,samplen,onset_samps);
        dataStruct(i).emg_epochs.(dataStruct(i).signals_f(j).name) = epochs;
        dataStruct(i).epoch_t = epoch_time;
    end
end
figure
rows = length(dataStruct);
cols = length(dataStruct(1).signals);
count = 1;
for i=1:rows
    epochData = dataStruct(i).emg_epochs;
    t = dataStruct(i).epoch_t;
    muscles = fieldnames(epochData);
    for j = 1:cols
        subplot(rows,cols, count);set(gca,'ButtonDownFcn',@fig_from_subplot);
        count =count+1;
        data = epochData.(muscles{j});
        hold on
        for p=1:size(data,1)
        plot(t,data(p,:));
        % stdev = std(dataStruct(i).signals_f(j).data);
        % ylim([-4*stdev, 4*stdev]);
        
        end
        hold off
        xlabel('time (s)');
        ylabel('amplitude (uV)')
        metadata = sprintf("%s numTrials %d, %s Stim, \n deltaF=%g Hz, amp=%s \n%s",dataStruct(i).signals(j).name,size(data,1),runMetaData(i).StimType, runMetaData(i).deltaF, runMetaData(i).Amplitude, runMetaData(i).Notes);
        title(metadata);
    end
end
figure
rows = length(dataStruct);
cols = length(dataStruct(1).signals);
count = 1;
for i=1:rows
    epochData = dataStruct(i).emg_epochs;
    t = dataStruct(i).epoch_t;
    muscles = fieldnames(epochData);
    for j = 1:cols
        subplot(rows,cols, count);set(gca,'ButtonDownFcn',@fig_from_subplot);
        count =count+1;
        data = epochData.(muscles{j});
        hold on
        avg = mean(data);
        
        stdev = std(data);
        c1 = avg+stdev;
        c2= avg-stdev;
        t2 = [t,fliplr(t)];
        inbetween = [c1,fliplr(c2)];
        fill(t2,inbetween,'r','EdgeColor','none','FaceAlpha',0.2)
        plot(t,avg,'r','LineWidth',3);
        % ylim([-4*stdev, 4*stdev]);
        
        
        hold off
        xlabel('time (s)');
        ylabel('amplitude (uV)')
        metadata = sprintf("%s numTrials %d, %s Stim, \n deltaF=%g Hz, amp=%s \n%s",dataStruct(i).signals(j).name,size(data,1),runMetaData(i).StimType, runMetaData(i).deltaF, runMetaData(i).Amplitude, runMetaData(i).Notes);
        title(metadata);
    end
end

%% Spectral Features
figure(5)
rows = length(dataStruct);
cols = length(dataStruct(1).signals);
count = 1;
for i=1:rows
    epochData = dataStruct(i).emg_epochs;
    muscles = fieldnames(epochData);
    for j = 1:cols
        subplot(rows,cols, count);set(gca,'ButtonDownFcn',@fig_from_subplot);
        count =count+1;
        data = epochData.(muscles{j});
        PXX = [];
        for p=1:size(data,1)
            [pxx,f] = pwelch((data(p,:)*10^(-6)),int32(dataStruct(i).fs/5),[],[],dataStruct(i).fs);
            % semilogy(f,pxx)
            PXX = [PXX pxx];
        end
        avg = mean(PXX')';
        
        stdev = std(PXX')';
        c1 = avg+stdev;
        c2= avg-stdev;
        t2 = [f,fliplr(f)];
        inbetween = [c1,fliplr(c2)];
        fill(t2,inbetween,'r','EdgeColor','none','FaceAlpha',0.2)
        loglog(f,avg,'r','LineWidth',3);
        hold on
    
        xlabel('time (s)');
        xlim([10 1000])
        ylabel('amplitude (uV)')
        metadata = sprintf("%s numTrials %d, %s Stim, \n deltaF=%g Hz, amp=%s \n%s",dataStruct(i).signals(j).name,size(data,1),runMetaData(i).StimType, runMetaData(i).deltaF, runMetaData(i).Amplitude, runMetaData(i).Notes);
        title(metadata);
    end
end

figure(6)
rows = length(dataStruct);
cols = length(dataStruct(1).signals);
count = 1;
for i=1:rows
    for j = 1:cols
        subplot(rows,cols, count);set(gca,'ButtonDownFcn',@fig_from_subplot);
        count =count+1;
        
        hold on
        data=dataStruct(i).signals_f(j).data;
        [pxx,f] = pwelch((data*10^(-6)),int32(dataStruct(i).fs),[],[],dataStruct(i).fs);
        hold off
        loglog(f,pxx)
        xlabel('time (s)');
        ylabel('amplitude (uV)')
        xlim([10 4100])
        metadata = sprintf("%s Run %d, %s Stim, \n deltaF=%g Hz, amp=%s \n%s",dataStruct(i).signals(j).name,runMetaData(i).Run,runMetaData(i).StimType, runMetaData(i).deltaF, runMetaData(i).Amplitude, runMetaData(i).Notes);
        title(metadata);
    end
end
%%
function [raw,filt] = process_EMG(input,channelLabels,fs)
raw = struct();
filt = struct();
sig1 = input(:,1:2:end);
sig2 = input(:,2:2:end);
for i=1:size(sig1,2)
    x = sig1(:,i)-sig2(:,i);
    raw(i).data = x - x(1);
    raw(i).name = channelLabels{i};
    filt(i).data = notch(bandpass(x-x(1),fs,20,450,3),fs,60,2);
    filt(i).name = channelLabels{i};
 
end
% for i=1:size(input,2):2
%     signals(i).data = input(:,i) - input(:,2*i+1);
%     signals(i).name = channelLabels{i};
% 
% end
end

function X=bandpass(data,fs,lowcut,highcut,order)
[b,a] = butter(order,[lowcut,highcut]/(fs/2));
X = filtfilt(b,a,data);
end

function X = reref_EEG(signals,chans)
X=0;
end

function locs = get_artifact_peaks(signal,peak_height,peak_distance,plotFlag)
 [pks,locs] = findpeaks(signal,"MinPeakHeight",peak_height,"MinPeakDistance",peak_distance);
     if plotFlag
         figure
         hold on
         plot(signal')
         scatter(locs,pks)
         hold off
     end
end

function X = epoch_by_peak_loc(data,peakLocs,samplen,onset_samples)
X=struct();
onset_samples = int64(onset_samples);
n_peaks = length(peakLocs);
trials = zeros(n_peaks,samplen);
for i=1:n_peaks
    if (samplen+peakLocs(i)-onset_samples) <= length(data) && peakLocs(i)-onset_samples+1 > 0
    trials(i,:) = data(peakLocs(i)-onset_samples+1:samplen+peakLocs(i)-onset_samples)';
    end
end
% X.trials=trials;
X=trials;
end

function X = notch(data,Fs,notchFreq,bandwidth)
% Parameters
% Sampling frequency in Hz
% Frequency to notch out in Hz
% Bandwidth around the notch frequency
d = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', ...
               notchFreq-bandwidth/2, 'HalfPowerFrequency2', ...
               notchFreq+bandwidth/2, 'DesignMethod', 'butter', 'SampleRate', Fs);
X = filtfilt(d,data);
end