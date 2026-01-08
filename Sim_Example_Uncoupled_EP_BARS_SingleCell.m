% Script to simulate 1 cell with the uncouple/recouple approach to
% producing an EP/BARS model in pseudo steady-state. 

clc
clear
close all

%% Load in necessary data
load('ISO100_X0BARS_SS','X0BARS');

%% Set parameters as needed 
param.verbose = 1; % give beat count
param.model = @model_Torord_dynCl_BARS_EP_Uncoupled_HF;
param.X0_Signaling = X0BARS; 
param.dt = 1; 
param.iso_dose = 1; % make sure this is congruent with ISO dose from X0BARS
param.bcl = 1000; % bcl in ms

%% Run simulation
tic
[time,X, currents] = run_ToR_ORd_BARS_EP_Uncoupled(param);
toc

%% Plot
V = [X{1}(:,1);X{2}(:,1);X{3}(:,1);X{4}(:,1);X{5}(:,1)];
t = (1:length(V))*param.dt;
figure;
plot(t,V);
xlabel('Time (ms)');
ylabel('V_m')
