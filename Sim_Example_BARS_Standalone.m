%% Script created to run the BARS standalone

clc
clear
close all

%% Run
param.iso_dose = 0.1; % dose of iso, muM
param.tspan = [0, 2e6]; % timespan, ms
[time,X] = run_BARS_standalone(param); % run the model

%% Plot state variables over entire timespan
for i = 1:size(X,2)
    hold on;
    plot(time,X(:,i),'-','color',[0 0 0 0.5]);
end

X0BARS = X(end,:); % save the end state

save('ISO010_X0BARS_SS','X0BARS');




