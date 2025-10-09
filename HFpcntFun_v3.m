function [HF, dCaRest,dCaiMax,dCaDelta, dCaD50,dAPD50,dAPD90, CaRestMat, CaiMaxMat, CaDeltaMat, CaD50Mat, APD50Mat, APD90Mat] = HFpcntFun_v3(bcl,WTbm,HFbm,pcnt)
 % calculate various biomarkers in one spot

% Add in the ratio of HF/WT:
kdiastolicCai = 165/96;
kCaimax = 367/746;
kCaiAmplitude = 398/804;
kCTD50 = 692/320;
%kCaTau = 306/209;
kAPD50 = 439.6/301.3;
kAPD90 = 304.3/233.2;

% Change between WT + HF
dCaRest = HFbm.diastolicCai/WTbm.diastolicCai;
dCaiMax = HFbm.Caimax/WTbm.Caimax;
dCaDelta = HFbm.CaiAmplitude/WTbm.CaiAmplitude;
dCaD50 = HFbm.CTD50/WTbm.CTD50;
dAPD50 = HFbm.APD50/WTbm.APD50;
dAPD90 = HFbm.APD90/WTbm.APD90;

% pre-alloc
CaRestMat = 0; 
CaiMaxMat = 0; 
CaDeltaMat = 0; 
CaD50Mat = 0; 
APD50Mat = 0; 
APD90Mat = 0; 


if dCaRest > (1+(kdiastolicCai-1)*pcnt/100); CaRestMat = 1; end
if dCaiMax < (1+(kCaimax-1)*pcnt/100); CaiMaxMat = 1; end
if dCaDelta < (1+(kCaiAmplitude-1)*pcnt/100); CaDeltaMat = 1; end
if dCaD50 > (1+(kCTD50-1)*pcnt/100); CaD50Mat = 1; end
if dAPD50 > (1+(kAPD50-1)*pcnt/100); APD50Mat = 1; end
if dAPD90 > (1+(kAPD90-1)*pcnt/100); APD90Mat = 1; end


if WTbm.APD90 < bcl & HFbm.DAD == 0 & HFbm.EAD == 0
    % Check HF condition
    if dCaRest > (1+(kdiastolicCai-1)*pcnt/100) && dCaiMax < (1+(kCaimax-1)*pcnt/100) && dCaDelta < (1+(kCaiAmplitude-1)*pcnt/100) && dCaD50 > (1+(kCTD50-1)*pcnt/100) && dAPD50 > (1+(kAPD50-1)*pcnt/100) && dAPD90 > (1+(kAPD90-1)*pcnt/100) % && dCaTau > (1+(kCaTau-1)*pcnt/100) 
        HF = 1;
    else
        HF = 0;
    end
else
    % Check HF condition, disregard diastolic Ca
    if dCaiMax < (1+(kCaimax-1)*pcnt/100) && dCaDelta < (1+(kCaiAmplitude-1)*pcnt/100) && dCaD50 > (1+(kCTD50-1)*pcnt/100) && dAPD50 > (1+(kAPD50-1)*pcnt/100) && dAPD90 > (1+(kAPD90-1)*pcnt/100) % && dCaTau > (1+(kCaTau-1)*pcnt/100) 
        HF = 1;
    else
        HF = 0;
    end
end

end