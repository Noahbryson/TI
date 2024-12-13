function modulation = modulationIndex(signal,envelope,windowlen)
% 
% 
% The depth of the amplitude modulation is computed by extracting the signals’
% envelope signals using Hilbert transform (abs(Hilbert), Matlab), and then dividing 
% half the peak-to-peak amplitude of the envelope signals by half the peak-to-peak 
% amplitude of the original signals.
% 
% 
% 
siglen = length(signal);
numFullSegments = floor(siglen/windowlen);
if mod(siglen,windowlen) == 0
    lastSegmentLength = windowlen;
else
    lastSegmentLength = windowlen + windowlen * (numFullSegments * windowlen < length(signal));
end
totalSegments = numFullSegments + (lastSegmentLength>windowlen);

% plot(signal),hold on, plot(envelope)
modulation = zeros(1,length(totalSegments));
for i=1:totalSegments
    startIdx = (i-1)*windowlen+1;
    endIdx = i*windowlen;
    if i == totalSegments
        endIdx = startIdx + lastSegmentLength -1;
    end
    % xline(startIdx,'r')
    % xline(endIdx,'b')
    env = envelope(startIdx:min(endIdx,length(envelope)));
    sig = signal(startIdx:min(endIdx,length(signal)));
    envPP = abs(max(env) - min(env))/2;
    sigPP = abs(max(sig) - min(sig))/2;

    modulation(i) = envPP/sigPP;
end

end