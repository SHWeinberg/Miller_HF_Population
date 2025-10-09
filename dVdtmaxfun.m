function [dVdtmax] = dVdtmaxfun(~,Vm,dt,~)

dV = diff(Vm);  % takes the maximum dV/dt
dVdt = dV./dt; % dVdt
dVdtmax = max(dVdt); % dVdt
dVdtmax = dVdtmax(1); 

end