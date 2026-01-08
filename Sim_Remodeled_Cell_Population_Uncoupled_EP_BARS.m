%% Script created to simulate a WT population using the ToR_ORd_BARS model

% This version has uncoupled EP & BARS activity. The model is recoupled
% halfway through the simulation. 

% Updating this script to generate and store data in individual cell files.
% This will result in a number of files equal to the size of the WT
% population.

clc
clear
close all

%% General Steps:
% 1. Load WT parameters
% 2. Generate HF parameters
% 3. Simulate HF cell
% 4. Save HF data
% 5. Go to next WT cell, generate HF pop, sim, and save

%% Section 1. Setting parameters
bcl = 1000;
ISO = 0.1; % be sure to change save name as appropriate at the bottomo of the script

%% Section 2. Load in datset
load('Control_Population_iso001'); paramsUncoupledCTRL = paramsUncoupledWT; % use correct naming convention
load('scaleIHF'); % use the same HF parameters as before.  
load('ISO010_X0SignalingSS'); 

%% Configure parallel cluster
% This section was removed because it was specific to the ohio
% supercomputer center. It is recommended that this process be
% parallelized via a cluster unless running small sample sizes.

%% Section 3. Generate HF multipliers (to be multiplied by each WT individual)
% Ranges were taken as HF ranges normalized to the WT mean

% Put indices/sizes
nCellstart = 1; % first cell simulated from the population of 1000
nCellend = 1e3; % last cell simulated from the population of 1000
nCell = nCellend - nCellstart + 1; % number of cells simulated
nHF = 1e3; % number of remodeled simulations

today = datetime('today'); 
dir = ''; % add directory save location here

% Pre-allocate
APD90 = nan(nHF,1);APD50 = nan(nHF,1);APD40 = nan(nHF,1);APD20 = nan(nHF,1);Tri9040 = nan(nHF,1);CTD50 = nan(nHF,1);CTD90 = nan(nHF,1);CaA = nan(nHF,1);CaiDiol = nan(nHF,1);RMP = nan(nHF,1);APA = nan(nHF,1);dVdtmax = nan(nHF,1);Vpeak = nan(nHF,1);EMwin = nan(nHF,1);Caimax = nan(nHF,1);avgCai = nan(nHF,1);X0NewHF = cell(nHF,1);
QNaCa = nan(nHF,1); QNa = nan(nHF,1);QNaL = nan(nHF,1);QKr = nan(nHF,1);QKs = nan(nHF,1);QK1 = nan(nHF,1);Qto = nan(nHF,1);QCaL = nan(nHF,1);QNaK = nan(nHF,1);QpCa = nan(nHF,1);QJup = nan(nHF,1);QNaCa_neg = nan(nHF,1);QNaCa_pos = nan(nHF,1);QpumpsE = nan(nHF,1);Nai = nan(nHF,1);Ki = nan(nHF,1);EAD = nan(nHF,1);DAD = nan(nHF,1);CaTSR = nan(nHF,1); CaNSR = nan(nHF,1); CaJSR = nan(nHF,1);
sus = nan(nHF,1); HF33 = nan(nHF,1); HF66 = nan(nHF,1); CaMKa = nan(nHF,1); 
fINap = nan(nHF,1); fICaLP = nan(nHF,1); fIKsP = nan(nHF,1); fPLBP = nan(nHF,1);
fTnIP = nan(nHF,1);fINaP = nan(nHF,1);fINaKP = nan(nHF,1);fRyRP = nan(nHF,1); fIKurP = nan(nHF,1);
cAMP_CAVV = nan(nHF,1);cAMP_ECAV = nan(nHF,1);cAMP_CYT = nan(nHF,1);inhib1_p = nan(nHF,1);

%% Section 4. Run these sims - loop through each cell
% Run through control dataset
for i = nCellstart:nCellend 

    tic 

    % Simulate the control to use as a comparison for HF change check
    [~,~, currentsControl] = run_ToR_ORd_BARS_EP_Uncoupled(paramsUncoupledCTRL(i));

    % start this individual's params
    paramsHF(1:nHF) = paramsUncoupledCTRL(i); % set scaling factors for control, next scale based on remodeling multipliers
      
    % Create the params structures
    for j = 1:nHF
        paramsHF(j).verbose = 1;
        paramsHF(j).iso_dose = ISO;
        paramsHF(j).maxTimePerBeat = 60; % max time per beat in sec.
        paramsHF(j).bcl = bcl;
        paramsHF(j).X0_Signaling = X0NewSignaling{j}; % set BARS signaling X0 based on HF remodeling case 
        
        % Multiply those that were scaled in WT
        paramsHF(j).INaL_Multiplier = paramsUncoupledCTRL(i).INaL_Multiplier*INaLHFMat(j);
        paramsHF(j).INaCa_Multiplier = paramsUncoupledCTRL(i).INaCa_Multiplier*INaCaHFMat(j);
        paramsHF(j).IK1_Multiplier = paramsUncoupledCTRL(i).IK1_Multiplier*IK1HFMat(j);
        paramsHF(j).INaK_Multiplier = paramsUncoupledCTRL(i).INaK_Multiplier*INaKHFMat(j);
        paramsHF(j).Jup_Multiplier = paramsUncoupledCTRL(i).Jup_Multiplier*JupHFMat(j);
        paramsHF(j).IpCa_Multiplier = paramsUncoupledCTRL(i).IpCa_Multiplier*IpCaHFMat(j);
        paramsHF(j).CaMKa_Multiplier =  paramsUncoupledCTRL(i).CaMKa_Multiplier*CaMKaHFMat(j);
        
        % Declare those that were not
        paramsHF(j).Jleak_Multiplier = paramsUncoupledCTRL(i).Jleak_Multiplier*JleakHFMat(j);
        paramsHF(j).tauhL_Multiplier = tauhLHFMat(j);
        paramsHF(j).Rb1_Multiplier = Rb1HFMat(j);
        paramsHF(j).ICaLdpHF_Multiplier = ICaLdpHFMat(j);
        paramsHF(j).PP1_Multiplier = PP1HFMat(j);
        paramsHF(j).PDE2_Multiplier = PDE2HFMat(j);
        paramsHF(j).PKA_Multiplier = PKAHFMat(j);
    end

    % Run simulation and save biomarkers
    for j = 1:nHF
        try % some simulations break, prevent code from being interrupted
        % Run Simulation
        [time,X, currentsHF] = run_ToR_ORd_BARS_EP_Uncoupled(paramsHF(j));
        % Save biomarkers
        [CaJSR(j), CaNSR(j), Nai(j), Ki(j), APD20(j), APD40(j), APD50(j), APD90(j), Tri9040(j), CTD50(j), CTD90(j), CaA(j), EMwin(j), Caimax(j), CaiDiol(j), avgCai(j), RMP(j), APA(j), dVdtmax(j), Vpeak(j), QNa(j), QNaL(j), QKr(j), QKs(j), QK1(j), Qto(j), QCaL(j), QNaK(j), QpCa(j), QJup(j), QNaCa(j), QNaCa_neg(j), QNaCa_pos(j), QpumpsE(j), DAD(j), EAD(j), X0NewHF{j},sus(j),HF33(j),HF66(j)] = funSaveBiomarkersHF_v2(paramsHF(j).dt,paramsHF(j).bcl,currentsHF,X,currentsControl);
        % Save phosphorylation - CaMKII
        CaMKa(j) = mean(currentsHF.CaMKa);
        fINap(j) = mean(currentsHF.fINap); % all p fractions are identical, can also be derived post-hoc via a simple equation using CaMKa
        % Save phosphorylation - PKA
        fICaLP(j) = mean(currentsHF.fICaLP);
        fIKsP(j) = mean(currentsHF.fIKsP);
        fPLBP(j) = mean(currentsHF.fPLBP);
        fTnIP(j) = mean(currentsHF.fTnIP);
        fINaP(j) = mean(currentsHF.fINaP);
        fINaKP(j) = mean(currentsHF.fINaKP);
        fRyRP(j) = mean(currentsHF.fRyRP);
        fIKurP(j) = mean(currentsHF.fIKurP);
        % Change other
        cAMP_CAVV(j) = mean([X{end-4}(1:end,69);X{end-3}(2:end,69);X{end-2}(2:end,69);X{end-1}(2:end,69);X{end}(2:end,69)]);
        cAMP_ECAV(j) = mean([X{end-4}(1:end,70);X{end-3}(2:end,70);X{end-2}(2:end,70);X{end-1}(2:end,70);X{end}(2:end,70)]);
        cAMP_CYT(j) = mean([X{end-4}(1:end,71);X{end-3}(2:end,71);X{end-2}(2:end,71);X{end-1}(2:end,71);X{end}(2:end,71)]);
        inhib1_p(j) = mean([X{end-4}(1:end,98);X{end-3}(2:end,98);X{end-2}(2:end,98);X{end-1}(2:end,98);X{end}(2:end,98)]);
        catch
           fprintf('\nCell Simulation Failed');
        end
    end
    
    timetaken = toc
    
    %% Subsection 4. Save data
    flnm = 'iso010_HF_Pop_ToRORdBARS_Epi_BCL'+string(paramsHF(1).bcl)+'_n1000_Cell'+string(i)+'_'+string(today);
    savenm = fullfile(dir,flnm);
    save(savenm,'X0NewHF','paramsHF','timetaken','APD90','APD40','APD50','APD20','CTD50','CTD90','APA','dVdtmax','CaA','CaiDiol','RMP','Tri9040','Vpeak','EMwin','Caimax','avgCai','CaJSR','CaNSR','QNa','QNaL','QKr','QKs','QK1','Qto','QCaL','QNaK','QpCa','QJup','QNaCa','QNaCa_neg','QNaCa_pos','QpumpsE','Nai','Ki','EAD','DAD','sus','HF33','HF66','CaMKa','fINap','fICaLP','fIKsP','fPLBP','fTnIP','fINaP','fINaKP','fRyRP','fIKurP','cAMP_CAVV','cAMP_ECAV','cAMP_CYT','inhib1_p'); 
    clear -vars paramsHF timetaken
end

