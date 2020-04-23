function [IF, IA, IP,CSC0]=InstSig(CSC,LFcutoff,HFcutoff,FiltOrder)

%2012-05-16 AndyP
%Bandpass filters signal then takes the Hilbert transform
% [IF, IA, IP, CSC0]=InstSig(CSC,LFcutoff,HFcutoff,FiltOrder);
% INPUTS
% CSC a tsd object
% LFcutoff [Hz] low frequency cutoff for bandpass filter
% HFcutoff [Hz] high frequency cutoff for bandpass filter
% OUTPUTS
% IF  tsd, range in [s], data in [Hz] instantaneous angular frequency
% IA  tsd, range in [s], data in [A.U.] instantaneous amplitude (power)
% IP  tsd, range in [s], data in [rad] instantaneous phase (radians)

assert(isa(CSC,'tsd'),'first input CSC must be a tsd object');
[CSC0]=FiltSig(CSC,LFcutoff,HFcutoff,FiltOrder);
[IF, IA, IP]=InstFreq(CSC0);
IF=tsd(CSC0.range,[IF(:); IF(end)]);
IA=tsd(CSC0.range,IA(:));
IP=tsd(CSC0.range,IP(:));
end




