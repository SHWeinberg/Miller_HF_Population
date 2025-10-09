% t = time 
% v = voltage or ion trace
% pcnt = percent repolarization
% win = 1/2 window to search upstroke
% thresh = detection threshold 
function [APD,tdown,downV,tup,upV] = APXfun(t,V,pcnt,thresh,win)

v1 = V(1:end-1); 
v2 = V(2:end); 
dt = t(2) - t(1);

% Phase 1, find upstrokes
upstrokei = find(v1 < thresh & v2 > thresh); % index voltage passes threshold (depolarization)
winind = cell(length(upstrokei),1); 
vwinup = cell(length(upstrokei),1); 
twinup = cell(length(upstrokei),1); 
vwindown = cell(length(upstrokei),1); 
twindown = cell(length(upstrokei),1); 
vmin = zeros(length(upstrokei),1); 
vmax = zeros(length(upstrokei),1); 
imin = zeros(length(upstrokei),1); 
imax = zeros(length(upstrokei),1); 
tmin = zeros(length(upstrokei),1); 
tmax = zeros(length(upstrokei),1); 
APA = zeros(length(upstrokei),1); 
upV = zeros(length(upstrokei),1); 
depoli = zeros(length(upstrokei),1); 
tup1 = zeros(length(upstrokei),1); 
tup2 = zeros(length(upstrokei),1); 
tup = zeros(length(upstrokei),1); 
vup1 = zeros(length(upstrokei),1); 
vup2 = zeros(length(upstrokei),1); 
downV = zeros(length(upstrokei),1); 
tdown1 = zeros(length(upstrokei),1); 
tdown2 = zeros(length(upstrokei),1); 
vdown1 = zeros(length(upstrokei),1); 
vdown2 = zeros(length(upstrokei),1); 
tdown = zeros(length(upstrokei),1); 
repoli = zeros(length(upstrokei),1); 
APD = zeros(length(upstrokei),1); 


% Phase 2, find exact time of upstroke 
for i = 1:length(upstrokei)
    if (upstrokei(i)-win/dt) > 0 && (upstrokei(i)+win/dt) < length(V)
        winind{i} = round((upstrokei(i)-win/dt):(upstrokei(i)+win/dt));
    elseif (upstrokei(i)+win/dt) > length(V)
        winind{i} = round((upstrokei(i)-win/dt):length(V));
    else
        winind{i} = round(1:(upstrokei(i)+win/dt));
    end

    vwinup{i} = V(winind{i}); % window of the upstroke
    twinup{i} = t(winind{i});
    [vmin(i), imin(i)] = min(vwinup{i}); % minimum of upstroke
    [vmax(i), imax(i)] = max(vwinup{i}); % maximum of upstroke
    tmin(i) = twinup{i}(imin(i)); 
    tmax(i) = twinup{i}(imax(i)); 

    APA(i) = vmax(i) - vmin(i); % action potential amplitude
    upV(i) = vmin(i) + 0.5*APA(i); % threshold for depolarization 
    depoli(i) = find(vwinup{i}>upV(i),1); % index just over threshold
    tup2(i) = twinup{i}(depoli(i));  % right over threshold
    tup1(i) = twinup{i}(max([depoli(i),2])-1); % right under threshold
    vup2(i) = vwinup{i}(depoli(i)); 
    vup1(i) = vwinup{i}(max([depoli(i),2])-1); 
    tup(i) = tup1(i) + (upV(i) - vup1(i))*(tup2(i) - tup1(i))/(vup2(i) - vup1(i)); %calculate interpolated upstroke time
end

% Phase 3, find repolarization times 
for i = 1:length(upstrokei)

    if i == length(upstrokei) %window to find repolarization 
        vwindown{i} = V(round(upstrokei(i)+win/dt:end)); % window is from this upstroke to the end of the AP
        twindown{i} = t(round(upstrokei(i)+win/dt:end)); % window is from this upstroke to the end of the AP
    else
        vwindown{i} = V(round(upstrokei(i)+win/dt: upstrokei(i+1)-win/dt)); % window is from this to next upstroke
        twindown{i} = t(round(upstrokei(i)+win/dt: upstrokei(i+1)-win/dt)); % window is from this to next upstroke
    end
    
    % repolarization threshold is the max voltage minus x% of APA
    downV(i) = vmax(i) - APA(i)*pcnt/100; % repolarization threshold is the max voltage minus x% of APA
    if sum(vwindown{i} < downV(i)) > 0
        repoli(i) = find(vwindown{i} < downV(i),1); 
    else
        APD(i) = nan; tup(i) = nan; tdown(i) = nan; upV(i) = nan; downV(i) = nan;
        break
    end
    tdown1(i) = twindown{i}(repoli(i));  % time before repolarizing
    tdown2(i) = twindown{i}(repoli(i)+1);  % time after repolarizing
    vdown1(i) = vwindown{i}(repoli(i)); % voltage before repoliarizing
    vdown2(i) = vwindown{i}(repoli(i)+1);  % voltage after repolarizing
    tdown(i) = tdown1(i) + (downV(i) - vdown1(i))*(tdown2(i) - tdown1(i))/(vdown2(i) - vdown1(i)); %calculate interpolated downstroke time
    
    APD(i) = tdown(i) - tup(i); 
end

if isempty(upstrokei)
    APD = nan; tdown = nan; tup = nan; downV = nan; upV  = nan;
end

if ~isempty(APD) & ~isnan(APD)
    APD = max(APD);
else
    APD = nan; 
end
    
    
    
    
    
    
    
    