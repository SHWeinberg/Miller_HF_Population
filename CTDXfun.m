function [CTDX] = CTDXfun(Time,Cat,dt,bcl,pcnt)
%function calculates CTDXX
    
 %loop over total number of beats
 for i = 1:5
    t_range = (1 + bcl*(i-1)*dt^-1):(bcl*i*dt^-1); %calculates time range for each beat
    t = Time(t_range); %pulls t values for time range
    cat = Cat(t_range); %pulls v values for time range
    
 % split AP into two using catmax 
    catmax = max(cat);
    catmax_ind = find(cat == catmax);
    
 % Calculate max and min 
    [catmax] = max(cat); %calculates the maximum cat
    [catmin] = min(cat); %calculates minimum cat
    catrange = catmax - catmin; %cat range (full APD cat range)
    CTDX = catmax - catrange*pcnt/100; %CTD X range
    CTD_50 = catmax - catrange*0.5; %CTD 50 range

 % This section used to find the upstroke time
    t1i = find(find(cat < CTD_50) < catmax_ind,1,'last'); %finds last number before -70 - index
    t2i = t1i + 1; %find first number after CTD50 - index
    cat1 = cat(t1i); %cat1
    cat2 = cat(t2i); %cat2
    t1 = t(t1i); %t1
    t2 = t(t2i); %t2
    cat3 = CTD_50; %upstroke cat
    
 if length(t1) < 1
     CTDX(i) = NaN;
 else
    t3 = t1 + (cat3 - cat1)*(t2 - t1)/(cat2 - cat1); %calculate interpolated time
    
 % This section used to calculate CTDXX - add component to calculate for
 % after max voltage
    t1i = find(cat > CTDX,1,'last'); %finds last number before CTDXX - index
    t2i = t1i + 1; %find first number after CTDXX - index
 if (t1i == length(cat)) || (t2i == length(cat))
    CTDX(i) = bcl;
 else
    cat1 = cat(t1i); %v1
    cat2 = cat(t2i); %v2
    t1 = t(t1i); %t1
    t2 = t(t2i); %t2
    cat4 = CTDX; %target voltage
    t4 = t1 + (cat4 - cat1)*(t2 - t1)/(cat2 - cat1); %calculate interpolated time
    
    CTDX(i) = t4 - t3; %calculate CTD by diff b/w upstroke and CTDXX (cat = CTDXX)
 end
 end
 end
 CTDX = max(CTDX); %max of CTDXX for the final 4 beats
 
 if CTDX > bcl
     CTDX = bcl;
 end
 
end