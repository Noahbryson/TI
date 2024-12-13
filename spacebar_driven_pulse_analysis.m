%% TI spacebar-pulse analysis
close all
user = expanduser('~'); % Get local path for interoperability on different machines, function in my tools dir.
if ispc
    boxpath = fullfile(user,"Box"); % Path to data
    BCI2KPath = "C:\BCI2000\BCI2000";
    addpath(genpath(fullfile(user,'Documents/NCAN/code/TI')))
else
    boxpath =  fullfile(user,'Library/CloudStorage/Box-Box'); % Path to data
    BCI2KPath = '/Users/nkb/Documents/NCAN/BCI2000tools';
    addpath(genpath(fullfile(user,'Documents/NCAN/code/TI')))
end

bci2ktools(BCI2KPath);
EEGkeys = {'FP1'; 'F3_233'; 'C3_234'; 'P3_235'; 'O1_236'; 'FP2'; 
    'F4_238';'C4_239';'P4_240';'O2_241';'F7_242';'T7_243';'P7';'F8_245'
    'T8_246';'P8';'F9_248';'F10_249';'FPZ';'FZ';'CZ';'PZ';'OZ'};
ECGkeys = {'EKG1';'EKG2'};
%% Loading and Epoching Data
datapath = fullfile(boxpath, 'Temporal-Interference-Project/data/');
% sessionPath = fullfile(datapath,'TI_LM1/20241112');
sessionPath = fullfile(datapath,'TI_SEEG_Tremor_BJH065/DAY6');
metadata = dir(fullfile(sessionPath,'_meta.csv'));
metadata = loadMetaData(fullfile(sessionPath,metadata.name));
tonicFlag = min(metadata{:,"duration_ms_"}) > 200; % tonic if the shortest duration is > 200ms
files = dir(fullfile(sessionPath,'*.dat'));
data = [];
rmChans = ismember(param.ChannelNames.Value,[EEGkeys;ECGkeys;{'REF1';'REF2'}]);
sig = sig(:,~rmChans);
param.ChannelNames.Value = param.ChannelNames.Value(~rmChans);
if tonicFlag
    for i=1:size(metadata,1)
        [sig,state,param] = load_bcidat(fullfile(sessionPath,metadata{i,'Fname'}{:}));
        info = struct;
        info.fs = param.SamplingRate.NumericValue;
        info.stimDuration = param.Stimuli.NumericValue(end,2);
        info.deltaF = metadata{i,'deltaF_Hz_'};
        info.baseF  = metadata{i,'base_Hz_'};
        info.amp    = metadata{i,'Amplitude_mA_'};
        temp = epochDatFile(sig,state,param,0,[],info);
        data = [data temp];
        
    end
else

    for i=1:size(metadata,1)
        [sig,state,param] = load_bcidat(fullfile(sessionPath,metadata{i,'Fname'}{:}));
        info = struct;
        info.fs = param.SamplingRate.NumericValue;
        info.stimDuration = param.Stimuli.NumericValue(end,2);
        info.deltaF = metadata{i,'deltaF_Hz_'};
        info.baseF  = metadata{i,'base_Hz_'};
        info.amp    = metadata{i,'Amplitude_mA_'};
        temp = epochDatFile(sig,state,param,0.2,[1 2;3 4;5 6]);
        data = [data temp];
    end
end



%% Avg Modulation Index Recruitment Curve 
avgMIs = aggregateTrialResultsAcrossChannels(data);
close all
figure(1)
ax = gca;
indat= {avgMIs.ModIdx};
legout=modulationRC(ax,[avgMIs.amp],[avgMIs.deltaF],indat);
legend(legout)

%% Modulation Index Recruitment Curve by Shank
figure(2)
shanks = unique({data.shank});
tcl = tiledlayout('flow');
for i=1:length(shanks)
avgMIs = aggregateTrialResultsAcrossChannels(data(ismember({data.shank},shanks{i})));
ax = nexttile(tcl);
indat= {avgMIs.ModIdx};
legout = modulationRC(ax,[avgMIs.amp],[avgMIs.deltaF],indat);
title(ax,strcat('Shank ',shanks{i}))
end
hL = legend(legout);
fontsize(hL,20,'points')
hL.Layout.Tile='East';



%%
function dataStruct = epochDatFile(data,states,params,postStimLen,diffLocs,info)
% postStimLen: post stimulation period added to epoch in seconds
data(isnan(data))=0;
chans = params.ChannelNames.Value;
fs = params.SamplingRate.NumericValue;


% get sample epochs
stimLen = str2num(regexprep(params.Stimuli.Value{end,2},'[a-zA-Z\s]', '')); %in seconds

sampleLen = floor((stimLen + postStimLen)*fs);
stimLocs = getInterval(states.StimulusCode, 2);
if postStimLen == 0
    postStimLocs = getInterval(states.StimulusCode, 3);
    tonic = 1;
end

% reshape and filter channels
dataStruct = struct;
if ~isempty(diffLocs) %if differential channel indexes are passed (rows: pairs, cols: channel idxs ie diffLocs(1,:) is one pair)

    for i=1:length(diffLocs) % handle the differential channels
        cname = regexprep(chans{diffLocs(i,:)}, '\d', '');
        dat = data(:,diffLocs(i,:));
        diffDat = dat(:,1) - dat(:,2);
        dataStruct(i).channel = cname;
        sig = butter_filter(diffDat,20,4,'high',fs);
        dataStruct(i).signal= sig;

        epochs = zeros(sampleLen,length(stimLocs));
        off = zeros(sampleLen,length(stimLocs));
        for j=1:length(stimLocs)
            epochs(:,j) = sig(stimLocs(j,1):stimLocs(j,1)+sampleLen-1);
            if tonic
                off(:,j) = sig(postStimLocs(j,1):postStimLocs(j,1)+sampleLen-1);
            end
        end
        dataStruct(i).epochs = epochs;
        if tonic
            dataStruct(i).offEpochs = off;
        end
        dataStruct(i).shank = cname;
    end
    % process single-ended channels
    rmLocs = unique(diffLocs(:));
    chans2 = chans(~rmLocs);
    data2 = data(:,~rmLocs);
    for i=1:length(chans2)
        sig = butter_filter(data2(:,i),1,2,'high',fs);
        dataStruct(i).signal= sig;
        dataStruct(i).channel = chans2{i};
        epochs = zeros(sampleLen,length(stimLocs));
        off = zeros(sampleLen,length(stimLocs));
        for j=1:length(stimLocs)
            epochs(:,j) = sig(stimLocs(j,1):stimLocs(j,1)+sampleLen-1);
            if tonic
                off(:,j) = sig(postStimLocs(j,1):postStimLocs(j,1)+sampleLen-1);
            end
        end
        dataStruct(i).epochs = epochs;
        if tonic
            dataStruct(i).offEpochs = off;
        end
        dataStruct(i).shank = regexprep(chans2{i},'\d','');
    end
else % handle only single ended channels ie just SEEG no EMG
    data = data - mean(data,2);
    for i=1:length(chans)
        cname = chans{i};
        locs = i;
        dat = data(:,locs);
        dataStruct(i).channel = cname;
        sig = butter_filter(dat,[1,200],2,'bandpass',fs);
        sig = notchFilt(sig,fs,1);
        env= abs(hilbert(sig));
        dataStruct(i).signal= sig;
        dataStruct(i).env= env;
        
        epochs = zeros(sampleLen,length(stimLocs));
        epochsEnv = zeros(sampleLen,length(stimLocs));
        onModIdxs = zeros(1,length(stimLocs));
        off = zeros(sampleLen,length(stimLocs));
        offEnv = zeros(sampleLen,length(stimLocs));
        offModIdxs = zeros(1,length(stimLocs));
        for j=1:length(stimLocs)
            epochs(:,j) = sig(stimLocs(j,1):stimLocs(j,1)+sampleLen-1);
            epochsEnv(:,j) = env(stimLocs(j,1):stimLocs(j,1)+sampleLen-1);
            onModIdxs(j) = mean(modulationIndex(epochs(:,j),epochsEnv(:,j),floor(5*fs/info.deltaF)));
            if tonic
                off(:,j) = sig(postStimLocs(j,1):postStimLocs(j,1)+sampleLen-1);
                offEnv(:,j) = env(postStimLocs(j,1):postStimLocs(j,1)+sampleLen-1);
                offModIdxs(j) = mean(modulationIndex(off(:,j),offEnv(:,j),floor(5*fs/info.deltaF)));
            end
        end
        dataStruct(i).epochs = epochs;
        dataStruct(i).epochsEnv = epochsEnv;
        dataStruct(i).ModIdx = onModIdxs;
        % dataStruct(i).ModIdxAv = mean(onModIdxs);
        if tonic
            dataStruct(i).offEpochs = off;
            dataStruct(i).offEpochsEn = offEnv;
            dataStruct(i).offModIdx = offModIdxs;
            % dataStruct(i).offModIdxAv = mean(offModIdxs);
        end
        dataStruct(i).shank = regexprep(chans{i},'\d','');
    end

end
[dataStruct(:).fs] = deal(info.fs);
[dataStruct(:).stimDuration] = deal(info.stimDuration); 
[dataStruct(:).deltaF] = deal(info.deltaF); 
[dataStruct(:).baseF] = deal(info.baseF);
[dataStruct(:).amp] = deal(info.amp);
end


function output = aggregateTrialResultsAcrossChannels(data)
amps = unique([data.amp]);
deltaF = unique([data.deltaF]);
[AA,BB] = meshgrid(amps,deltaF);
combos = [AA(:) BB(:)];
if isfield(data,'offModIdx')
    offModFlag = 1;
else
    offModFlag =0;
end
output = struct;
for i=1:size(combos,1)
ampLoc = ismember([data.amp],combos(i,1));
fLoc = ismember([data.deltaF],combos(i,2));
selLoc = ampLoc & fLoc; % where target values are true
entries = data(selLoc);
ModIdx = mean(cell2mat({entries.ModIdx}'),1);
output(i).amp = combos(i,1);
output(i).deltaF = combos(i,2);
output(i).ModIdx = ModIdx;
if offModFlag
    offModIdx = mean(cell2mat({entries.offModIdx}'),1);
    output(i).offModIdx = offModIdx;
end

end
end


function legOut = modulationRC(ax,ampArray,fArray,MI)
% ylim(ax,[0 2]);
deltaFs = unique(fArray);
colors = {'r' 'b', 'g','p'};
hold on
data = zeros(length(MI),1);
stdev = zeros(length(MI),1);
for i=1:length(MI)
    data(i) = mean(MI{i});
    stdev(i) = std(MI{i});
end
legOut = cell(length(deltaFs)*2,1);
for i=1:length(deltaFs)
loc = ismember(fArray,deltaFs(i));
x = ampArray(loc);
y = data(loc);
errorbar(ax,x,y,stdev(loc),'LineStyle','none','Color',colors{i},'DisplayName','')
scatter(ax,x,y,'Color',colors{i},'DisplayName',sprintf('delta F: %0.1f Hz',deltaFs(i)));
legOut{2*i-1} = '';
legOut{2*i} = sprintf('delta f: %0.1f', deltaFs(i));
% leg = arrayfun(@num2str,deltaFs,'UniformOutput',false);

end
% legend show
end

function metadata = loadMetaData(fp)
metadata = readtable(fp);
end


function intervals = getInterval(stimCode,value)
A = find(stimCode==value);
differences = diff(A);

% Find indices where the difference is greater than 1
if ~isempty(differences)
    breaks = find(differences > 1);
    intervals = zeros(length(breaks),2);
    % Define start and end indices of each interval
    idx = A(1);
    for i=1:length(breaks)
        intervals(i,:) = [idx,A(breaks(i))];
        idx = A(breaks(i)+1);
    end
    % disp(intervals);
else
    intervals = [];
end
end