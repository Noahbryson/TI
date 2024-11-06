clear 
user = expanduser('~'); % Get local path for interoperability on different machines, function in my tools dir. 
if ispc
    boxpath = fullfile(user,'Box/'); % boxpath
    BCI2KPath = 'C:\BCI2000\BCI2000';
else
    boxpath =  fullfile(user,'Library/CloudStorage/Box-Box'); % boxpath
    BCI2KPath = '/Users/nkb/Documents/NCAN/BCI2000tools';
end

savepath = fullfile(boxpath,'Temporal-Interference-Project/BCI2000/Keysight');
bci2ktools(BCI2KPath);
%%
deviceAddress = 'USB0::10893::36097::CN62200145::0::INSTR'; 
% deviceAddress: USB VISA Address for our keysight,
%       address must be changed for other hardware but 
%       seems static across computers for our device.
param = struct;
functionComment = 'Shape of ch1 waveform: 1 Sinusoid, 2 Square, 3 Triangle, 4 Ramp, 5 Pulse, 6 Noise, 7 DC, 8 PRBS, 9 User (enumeration)';
Voltage = 1; % V
DC_offset = 0; % V
baseFreq = 10000; % hz
deltaF = 20; % hz
trigger1 = 37; % left arrow 
trigger2 = 39; % right arrow
triggerOff = 40; % down arrow
VoltageLims = 0; %V, set to zero
triggerState = 'KeyUp';
stimExperimentFlag = 1;
numTrials = 10;
if stimExperimentFlag
    triggerState = 'StimulusCode';
    trigger1 = 2; % stimulation on stim code
    trigger2 = 2; % stimulation on stim code
    triggerOff = 3; % stimulation off stim code
    
end
% % % % % ch1
param = writeFields(param,'Address','Keysight','string','','','','',deviceAddress);
param = writeFields(param,'UsingCh1','Keysight','int','','','','',1);
param = writeFields(param,'Frequency1','Keysight','float','','','','',baseFreq);
param = writeFields(param,'Voltage1','Keysight','float','','','','Signal Amplitude',Voltage);
param = writeFields(param,'VoltageOffset1','Keysight','float','','','','DC Offset of Signal',DC_offset);
param = writeFields(param,'HighVoltage1','Keysight','float','','','','',VoltageLims);
param = writeFields(param,'LowVoltage1','Keysight','float','','','','',-VoltageLims);
param = writeFields(param,'TriggerOnState1','Keysight','string','','','','',triggerState);
param = writeFields(param,'TriggerOnValue1','Keysight','int','','','','',trigger1);
param = writeFields(param,'TriggerOffState1','Keysight','string','','','','',triggerState);
param = writeFields(param,'TriggerOffValue1','Keysight','int','','','','',triggerOff);
param = writeFields(param,'Phase1','Keysight','float','0','-360','360','Signal Phase (-360 to 360)',0);
param = writeFields(param,'Load1','Keysight','float','50Ohm','','','','50Ohm');
param = writeFields(param,'Function1','Keysight','int','','','',functionComment,1); 
% % % % % ch2
param = writeFields(param,'UsingCh2','Keysight','int','','','','',1);
param = writeFields(param,'Frequency2','Keysight','float','','','','',baseFreq+deltaF);
param = writeFields(param,'Voltage2','Keysight','float','','','','Signal Amplitude',Voltage);
param = writeFields(param,'VoltageOffset2','Keysight','float','','','','DC Offset of Signal',DC_offset);
param = writeFields(param,'HighVoltage2','Keysight','float','','','','',VoltageLims);
param = writeFields(param,'LowVoltage2','Keysight','float','','','','',-VoltageLims);
param = writeFields(param,'TriggerOnState2','Keysight','string','','','','',triggerState);
param = writeFields(param,'TriggerOnValue2','Keysight','int','','','','',trigger2);
param = writeFields(param,'TriggerOffState2','Keysight','string','','','','',triggerState);
param = writeFields(param,'TriggerOffValue2','Keysight','int','','','','',triggerOff);
param = writeFields(param,'Phase2','Keysight','float','0','-360','360','Signal Phase (-360 to 360)',0);
param = writeFields(param,'Load2','Keysight','float','50Ohm','','','','50Ohm');
param = writeFields(param,'Function2','Keysight','int','','','',functionComment,1); 

if stimExperimentFlag

    stimRowLabs = {'caption';'icon';'av';'EarlyOffsetExpression';'StimulusDuration'};
    param.Stimuli.RowLabels    = stimRowLabs;
    param.Stimuli.Section      = 'Application';
    param.Stimuli.Type         = 'matrix';
    param.Stimuli.DefaultValue = '';
    param.Stimuli.LowRange     = '';
    param.Stimuli.HighRange    = '';
    param.Stimuli.Comment      = 'captions and icons to be displayed, sounds to be played for different stimuli';
    param.Stimuli.Value        = cell(5,4); % 

    param.Stimuli.Value{1,1}   = 'Press Space To Stimulate'; 
    param.Stimuli.Value{2,1}   = '';
    param.Stimuli.Value{3,1}   = '';
    param.Stimuli.Value{4,1}   = 'Keydown==32';
    param.Stimuli.Value{5,1}   = '40s';

    param.Stimuli.Value{1,trigger1}   = ''; 
    param.Stimuli.Value{2,trigger1}   = '';
    param.Stimuli.Value{3,trigger1}   = '';
    param.Stimuli.Value{4,trigger1}   = '';
    param.Stimuli.Value{5,trigger1}   = '0.1s';

    param.Stimuli.Value{1,triggerOff}   = 'Stim Off, Please Wait'; 
    param.Stimuli.Value{2,triggerOff}   = '';
    param.Stimuli.Value{3,triggerOff}   = '';
    param.Stimuli.Value{4,triggerOff}   = '';
    param.Stimuli.Value{5,triggerOff}   = '5s';

    param.Stimuli.Value{1,4}   = 'End of Run'; 
    param.Stimuli.Value{2,4}   = '';
    param.Stimuli.Value{3,4}   = '';
    param.Stimuli.Value{4,4}   = '';
    param.Stimuli.Value{5,4}   = '1s';
    % sequence type
    param.SequenceType.Section       = 'Application';
    param.SequenceType.Type          = 'int';
    param.SequenceType.DefaultValue  = '';
    param.SequenceType.LowRange      = '';
    param.SequenceType.HighRange     = '';
    param.SequenceType.Comment       = 'Sequence of stimuli is 0 deterministic, 1 random, 2 P3Speller compatible (enumeration)';
    param.SequenceType.Value         = {'0'};
    % sequence
    param.Sequence.Section       = 'Application';
    param.Sequence.Type          = 'intlist';
    param.Sequence.DefaultValue  = '';
    param.Sequence.LowRange      = '';
    param.Sequence.HighRange     = '';
    param.Sequence.Comment       = '';
    param.Sequence.Value         = cell(numTrials*3+1,1);
    c = 1;
    for i=1:numTrials
        param.Sequence.Value{c,1} = '1';
        c = c+1;
        param.Sequence.Value{c,1} = '2';
        c = c+1;
        param.Sequence.Value{c,1} = '3';
        c = c+1;
    end
    param.Sequence.Value{c,1} = '4';
end






% Write 2 file
fname = sprintf('%.1fV_%.1fhz_%.1fhz-deltaf_keysightTI.prm',Voltage,baseFreq,deltaF);
if stimExperimentFlag
    fname = strcat('StimExp_',fname);
end
filename = fullfile(savepath,fname);
parameter_lines = convert_bciprm( param );
fid = fopen(filename, 'w');

for loc=1:length(parameter_lines)
    fprintf( fid, '%s', parameter_lines{loc} );
    fprintf( fid, '\r\n' );
end
fclose(fid);
%%
function param = writeFields(param,fieldname,section,Type,default,low,high,comment,value)
% function to quickly write single values to param file. will not work in
% current implementation for matrix type params
% if isa(value,'int') | isa(value,'double')
% 
% end
if isnumeric(value)
    value = num2str(value);
end
value = {value};
param.(fieldname).Section = section;
param.(fieldname).Type = Type;
param.(fieldname).DefaultValue = default;
param.(fieldname).LowRange = low;
param.(fieldname).HighRange = high;
param.(fieldname).Comment = comment;
param.(fieldname).Value = value;

end
