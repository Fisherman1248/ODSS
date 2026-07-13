% ODSS code, generated from offical OTFS matlab example.
clc; clear; close all;

%% Simulation Set up
q = 2;
Nscale = 7;             % n = 0,...,6, as in the paper example
B = 1280;               % system bandwidth in Hz, fig. 8 illustrate
W = B*(q-1)/(q^Nscale-1); % Transmit filter bandwidth, Eq. (77)
T = 1.9;                % ODSS pulse/block duration in seconds


%% Grid Population for Channel Sounding
% Pilot generation and grid population
pilotBin = floor(N/2)+1;
Pdd = zeros(M,N);
Pdd(1,pilotBin) = exp(1i*pi/4); % populate just one bin to see the effect through the channel

%% OTFS Modulation
txOut = helperOTFSmod(Pdd,padLen,padType);
