function param = keysightParmGen(V,DC,func1,func2,freq1,freq2,stimState,stimOn,StimOff)
param = struct;

param.UsingCh{1} = 1;
param.Load{1} = 50;
param.Function{1} = func1;
param.Phase{1} = 0;
param.Frequency{1} = freq1;
param.Voltage{1} = V;
param.VoltageOffset{1} = DC;
param.HighVoltage{1} = 0;
param.LowVoltage{1} = 0;
param.TriggerOnState{1} = stimState;
param.TriggerOnValue{1} = stimOn;
param.TriggerOffState{1} = StimState;
param.TriggerOffValue{1} = stimOff;



end