function [BIN_EAD] = EADfun(Vm,dt,bcl,APD90)
% Check if an EAD occurs by seeing if there is a local minimum between the
% initial notch and the APD90.

if ~isnan(APD90)
    % Preallocate
    local_mins = zeros(1,5);

    time = (1:length(Vm))*dt';

    % loop through each beat, checking for EADs
    for i = 1:5
        t_range = (1 + bcl*(i-1)*dt^-1):(bcl*i*dt^-1); %calculates time range for each beat
        t = time(t_range); %pulls t values for time range
        v = Vm(t_range); %pulls v values for time range

        t1i = 50/dt; %find(t == 20); %t1 is at 20 ms
        t2i = round(APD90)/dt; %find(t == round(APD90,1)+4.5e4+1e3*(i-1)); %time of APD90
        local_mins(i) = sum(islocalmin(v(t1i:t2i)) == 1); % find local min
    end

    if sum(local_mins) > 0
        BIN_EAD = 1;
    else
        BIN_EAD = 0;
    end
else
    BIN_EAD = 1;
end


end