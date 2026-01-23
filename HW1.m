% Q5

clear; close all; clc;

N = 100;
Ns = 10; % Number of symbols to be plotted
Tsym = 1; % Symbol time in seconds
oversampleFactor = 100;
Tres = Tsym / oversampleFactor;
t = 0:Tres:Ns*Tsym; % Time axis for symbols to be plotted

bc = sign(2 * rand(1,N) - 1);
bs = sign(2 * rand(1,N) - 1);



% uc = zeros(1, size(t));
% us = zeros(1, size(t));
% 
% for i = 1:Ns
% 
% 
% end

