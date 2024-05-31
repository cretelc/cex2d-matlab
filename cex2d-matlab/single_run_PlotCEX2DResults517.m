% Single run PlotCEX2DResults517.m on all Output files 
close all; clear all; clc;

% Prompt user to identify the file to assess.
filename = uigetfile('*.txt','Select CEX2D output file');
% Display the chosen file. 
disp(filename);
% Run the plotting code. 
cex2dresults_function(filename);

%% Functions 