%% subgroup discovery via Quadratic Programming
% reference:
% Li, R, Perneczky, R, Drzezga, A, and Kramer, S (2013).
% Efficient Redundancy Reduced Subgroup Discovery via Quadratic Programming
% accepted; Journal of Intelligent Information Systems, 
% invited as best Discovery Science 2012 papers.

clc
clear
close all

load pima

thres = 0.01;
beamWidth = 5;
alpha = 0.5;

%% discretize the continuous data
[data cuts] = entroDis(data, label);

%% quadratic programming approach
SDrules1 = sdQPMIFS(data, label, thres, alpha);

%% beam search
SDrules2 = sdBeam(data, label, thres, beamWidth);

%% optimistic estimate
SDrules3 = tOptiEst(data, label, thres);


