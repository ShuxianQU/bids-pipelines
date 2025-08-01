function [b0map, freqmap] = calculate_b0map (phamap, dte)
% 
% CALCULATE_B0MAP Calculate B0 maps based on a phase map in Siemens format. 
%
% Usage: [b0map, freqmap] = calculate_b0map (phamap, dte)
%
% Returns
% -------
% 
% b0map: b0 map in uT.
% freqmap: b0 map in Hz.
% 
%
% Expects
% -------
% phamap: phase map in Siemens magic numbers. 
% dte: TE difference (in ms) between two echos.
%
%
% See also: calculate_b0map_dicom
%
%
% Copyright (C) 2020 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Thu Jan 30 16:29:57 2020
%

dGamma      = 42576375;
dT= dte* 1e-3; % sec

% convert
b0map = 1e6* (double(phamap)-2048.0)/4096.0/dGamma/dT;
freqmap = (double(phamap)-2048.0)/4096.0/dT;

disp('-> Done...')
