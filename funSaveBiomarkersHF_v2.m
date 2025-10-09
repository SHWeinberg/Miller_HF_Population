
function [CaJSR, CaNSR, Nai, Ki, APD20, APD40, APD50, APD90, Tri9040, CTD50, CTD90, CaiAmplitude, EMwin, Caimax, diastolicCai, avgCai, RMP, APA, dVdtmax, Vpeak, QNa, QNaL, QKr, QKs, QK1, Qto, QCaL, QNaK, QpCa, QJup, QNaCa, QNaCa_neg, QNaCa_pos, QpumpsE, DAD, EAD, X0NewHF,sus,HF33,HF66, dCaRest,dCaiMax,dCaDelta, dCaD50,dAPD50,dAPD90] = funSaveBiomarkersHF_v2(dt,bcl,currentsHF,X,currentsWT)
% Edited from v1 to insert WT currents to calculate biomarkers for
% comparison when calculating WT->HF changes
try
    if ~isempty(X)
        % Calculate concentrations
        CaJSR = max(currentsHF.CaJSR);
        CaNSR = max(currentsHF.CaNSR);
        Nai = max([X{1}(:,2), X{2}(:,2), X{3}(:,2), X{4}(:,2), X{5}(:,2)],[],'all');
        Ki = max([X{1}(:,4), X{2}(:,4), X{3}(:,4), X{4}(:,4), X{5}(:,4)],[],'all');

        % Calculate biomarkers
        APD20 = min([1000,APXfun(currentsHF.time,currentsHF.V,20,-50,5)]);
        APD40 = min([1000,APXfun(currentsHF.time,currentsHF.V,40,-50,10)]);
        APD50 = min([1000,APXfun(currentsHF.time,currentsHF.V,50,-50,10)]);
        APD90 = min([1000,APXfun(currentsHF.time,currentsHF.V,90,-50,10)]);
        Tri9040 = APD90-APD40;
        CTD50 = min([1000,CTDXfun(currentsHF.time,currentsHF.Cai,dt,bcl,50)]);
        CTD90 = min([1000,CTDXfun(currentsHF.time,currentsHF.Cai,dt,bcl,90)]);
        [CaiAmplitude, diastolicCai, Caimax]= CaiAmplitudefun(currentsHF.Cai,dt,bcl);
        %CaTau = CaTauFun(currentsHF.Cai,dt,bcl);
        EMwin = CTD90 - APD90; % electromechanical window
        avgCai = mean(currentsHF.Cai);
        RMP = min(currentsHF.V);
        APA = APAfun(currentsHF.time,currentsHF.V,dt,bcl);
        dVdtmax = dVdtmaxfun(currentsHF.time,currentsHF.V,dt,bcl);
        Vpeak = max(currentsHF.V);

        % Calculate current integrals
        QNa = abs(trapz(currentsHF.time,currentsHF.INa));
        QNaL = abs(trapz(currentsHF.time,currentsHF.INaL));
        QKr = abs(trapz(currentsHF.time,currentsHF.IKr));
        QKs = abs(trapz(currentsHF.time,currentsHF.IKs));
        QK1 = abs(trapz(currentsHF.time,currentsHF.IK1));
        Qto = abs(trapz(currentsHF.time,currentsHF.Ito));
        QCaL = abs(trapz(currentsHF.time,currentsHF.ICaL));

        % Calc Pump/Exchanger integral
        QNaK = abs(trapz(currentsHF.time,currentsHF.INaK));
        QpCa = abs(trapz(currentsHF.time,currentsHF.IpCa));
        QJup = abs(trapz(currentsHF.time,currentsHF.Jup));

        % Positive vs. negative NCX
        INaCa_neg = currentsHF.INaCa;
        INaCa_neg(INaCa_neg > 0) = 0; % remove positive component
        INaCa_pos = currentsHF.INaCa;
        INaCa_pos(INaCa_pos < 0) = 0; % remove negative component
        QNaCa = abs(trapz(currentsHF.time,currentsHF.INaCa));
        QNaCa_neg = abs(trapz(currentsHF.time,INaCa_neg));
        QNaCa_pos = abs(trapz(currentsHF.time,INaCa_pos));

        % This needs to be done to consider the charge/ATP
        % INaK = +1 charge per 1ATP
        % IpCA = +2 charge per 1ATP (straight up based on the wikipedia page)
        % SERCA = +4 charge per 1ATP (https://doi.org/10.1186%2Fs13395-021-00280-7)
        QpumpsE = abs(QNaK) + abs(QpCa)/2 + abs(QJup)/4; % total pump activity

        DAD = DADfun(currentsHF.V,dt,bcl,APD90);
        EAD = EADfun(currentsHF.V,dt,bcl,APD90);

        % feed into the HF phenotype check function
        HFbm.DAD = DAD;
        HFbm.EAD = EAD; 

        %% Add in susceptibility or HF case check

        sus = SusCheckFun(APD90,DAD,EAD);

        % biomarkers needed: Ca Amplitude, diastolic Ca, Ca max, APD90(?),
        % CTD50, Ca Tau (needs to be calculated)
        % Create structures to input to HF function

        % WT
        [WTbm.CaiAmplitude, WTbm.diastolicCai, WTbm.Caimax]= CaiAmplitudefun(currentsWT.Cai,dt,bcl);
        WTbm.CTD50 = min([1000,CTDXfun(currentsWT.time,currentsWT.Cai,dt,bcl,50)]);
        %WTbm.CaTau = CaTauFun(currentsWT.Cai,dt,bcl);
        WTbm.APD90 = min([1000,APXfun(currentsWT.time,currentsWT.V,90,-50,10)]);
        WTbm.APD50 = min([1000,APXfun(currentsWT.time,currentsWT.V,50,-50,10)]);

        % HF
        HFbm.CaiAmplitude = CaiAmplitude;
        HFbm.diastolicCai = diastolicCai;
        HFbm.Caimax = Caimax;
        HFbm.CTD50 = CTD50;
        HFbm.APD50 = APD50;
        HFbm.APD90 = APD90;

        % Calculate binary classification

        [HF33] = HFpcntFun_v3(bcl,WTbm,HFbm,33);
        [HF66, dCaRest,dCaiMax,dCaDelta, dCaD50,dAPD50,dAPD90] = HFpcntFun_v3(bcl,WTbm,HFbm,66);

    else
        % Calculate concentrations
        CaJSR = nan;
        CaNSR = nan;
        Nai = nan;
        Ki =nan;

        % Calculate biomarkers
        APD20 = nan;
        APD40 =nan;
        APD50 = nan;
        APD90 = nan;
        Tri9040 = nan;
        CTD50 = nan;
        CTD90 = nan;
        CaiAmplitude = nan;
        EMwin = nan;
        Caimax = nan;
        diastolicCai = nan;
        avgCai = nan;
        RMP = nan;
        APA =nan;
        dVdtmax =nan;
        Vpeak = nan;

        % Calculate current integrals
        QNa = nan;
        QNaL = nan;
        QKr = nan;
        QKs = nan;
        QK1 = nan;
        Qto =nan;
        QCaL = nan;

        % Calc Pump/Exchanger integral
        QNaK = nan;
        QpCa = nan;
        QJup = nan;

        % Positive vs. negative NCX
        QNaCa = nan;
        QNaCa_neg = nan;
        QNaCa_pos = nan;

        % This needs to be redone to consider the charge/ATP
        % INaK = +1 charge per 1ATP
        % IpCA = +2 charge per 1ATP (straight up based on the wikipedia page)
        % SERCA = +4 charge per 1ATP (https://doi.org/10.1186%2Fs13395-021-00280-7)
        QpumpsE = nan;

        DAD = nan;
        EAD = nan;

        sus = nan;
        HF33 = nan;
        HF66 = nan;

        dCaRest = nan;
        dCaiMax = nan;
        dCaDelta = nan;
        dCaD50 = nan;
        dAPD50 = nan;
        dAPD90 = nan;


    end

    % Handle cells that were terminated prematurely
    if ~isempty(X)
        % Save end-state vectors
        X0NewHF = X{end-4}(1,:); % save the end state for each trace as a vector. Then I can easily simulate any additional data
    else
        X0NewHF = nan;
    end

catch
    fprintf('Error calculating biomarkers, replacing with nan');
    % Calculate concentrations
    CaJSR = nan;
    CaNSR = nan;
    Nai = nan;
    Ki =nan;

    % Calculate biomarkers
    APD20 = nan;
    APD40 =nan;
    APD50 = nan;
    APD90 = nan;
    Tri9040 = nan;
    CTD50 = nan;
    CTD90 = nan;
    CaiAmplitude = nan;
    EMwin = nan;
    Caimax = nan;
    diastolicCai = nan;
    avgCai = nan;
    RMP = nan;
    APA =nan;
    dVdtmax =nan;
    Vpeak = nan;

    % Calculate current integrals
    QNa = nan;
    QNaL = nan;
    QKr = nan;
    QKs = nan;
    QK1 = nan;
    Qto =nan;
    QCaL = nan;

    % Calc Pump/Exchanger integral
    QNaK = nan;
    QpCa = nan;
    QJup = nan;

    % Positive vs. negative NCX
    QNaCa = nan;
    QNaCa_neg = nan;
    QNaCa_pos = nan;

    % This needs to be redone to consider the charge/ATP
    % INaK = +1 charge per 1ATP
    % IpCA = +2 charge per 1ATP (straight up based on the wikipedia page)
    % SERCA = +4 charge per 1ATP (https://doi.org/10.1186%2Fs13395-021-00280-7)
    QpumpsE = nan;

    DAD = nan;
    EAD = nan;
    X0NewHF = nan;

    sus = nan;
    HF33 = nan;
    HF66 = nan;

    dCaRest = nan;
    dCaiMax = nan;
    dCaDelta = nan;
    dCaD50 = nan;
    dAPD50 = nan;
    dAPD90 = nan;

end

