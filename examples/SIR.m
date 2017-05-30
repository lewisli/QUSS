close all; clear all;
addpath('../src/evidential/');

load('../data/SIR/3dslresults/Prior.mat');
load('../data/SIR/3dslresults/Time1000/Observed.mat')

FontSize = 30;
PlotResponses(PriorData,TrueData,FontSize);

load('../data/SIR/parameters/Prior/Sampled.mat');
