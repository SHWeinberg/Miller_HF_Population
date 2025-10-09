function [BIN_DAD] = DADfun(Vm,dt,bcl,APD90)

% Preallocate
local_maxes = zeros(1,5);
time = (1:length(Vm))*dt';

% loop through each beat, checking for DADs
 for i = 1:5
    t_range = round(1 + bcl*(i-1)*dt^-1):round(bcl*i*dt^-1); %calculates time range for each beat
    t = time(t_range); %pulls t values for time range
    v = Vm(t_range); %pulls v values for time range

   t1i = round(APD90)/dt; %find(t == round(APD90,1)+4.5e4+1e3*(i-1)); %time of APD90
   t2i = round(bcl)/dt; %find(t == 20); %t1 is at 20 ms
   local_maxes(i) = sum(islocalmax(v(t1i:t2i),'MinProminence',5) == 1); % find local min
 end

 if sum(local_maxes) > 0 % greater than 0
     BIN_DAD = 1;
 else
     BIN_DAD = 0;
 end

end
