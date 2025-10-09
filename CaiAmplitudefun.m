function [CaAmplitude, CaDiol, CaMax] = CaiAmplitudefun(Cat,dt,bcl)
%function calculates Cai amplitude
CaAmplitude = zeros(1,5);
%loop over total number of beats
for i = 1:5
    t_range = (1 + bcl*(i-1)*dt^-1):(bcl*i*dt^-1); %calculates time range for each beat
    cat = Cat(t_range); %pulls v values for time range

    % Calculate max and min
    [catmax(i)] = max(cat); %calculates the maximum cat
    [catmin(i)] = min(cat); %calculates minimum cat
    CaAmplitude(i) = catmax(i) - catmin(i); %cat range (full APD cat range)
end
CaDiol = max(catmin);
CaMax = min(catmax);
CaAmplitude = max(CaAmplitude(i)); %max of CTD90 for the final 4 beats
end