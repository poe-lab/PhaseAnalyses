function [IF_theta, IA_theta, IP_theta, CSC0_theta]=InstSig_theta(CSC, varargin)
%2012-05-16 AndyP

%Bandpass filters signal then takes the Hilbert transform
% [IF, IA, IP, CSC0]=InstSig_theta(CSC);
% 
% calls InstSig with 
%   LFcutoff = 6
%   HFcutoff = 10;
%   FiltOrder = 256;

LFcutoff = 6; %[Hz]
HFcutoff = 10; %[Hz]
FiltOrder = 256;
process_varargin(varargin);
[IF_theta , IA_theta , IP_theta, CSC0_theta ]=InstSig(CSC,LFcutoff,HFcutoff,FiltOrder);
end




