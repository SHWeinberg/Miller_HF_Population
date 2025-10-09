function [APA] = APAfun(~,Vm,dt,bcl)
%function calculates max Voltage
 %loop over total number of beats
 
 Vmax = zeros(1,5);
 rmp = zeros(1,5);
 APA = zeros(1,5);
 
 for i = 1:5
    t_range = ceil((1 + bcl*(i-1)*dt^-1)):ceil((bcl*i*dt^-1));
    v = Vm(t_range);
    Vmax(i) = max(v); 
    rmp(i) = min(v);
    APA(i) = abs(Vmax(i) - rmp(i)); 
 end
 APA = max(APA);
end