function [sus] = SusCheckFun(APD90,DAD,EAD)
% Check if the cell has an EAD, DAD, APD90 > 500

%if APD90 > 500 || DAD == 1 || EAD == 1
if APD90 > 500 || EAD == 1
    sus = 1;
else 
    sus = 0;
end


end