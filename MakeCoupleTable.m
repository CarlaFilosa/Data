function [TLabels] = MakeCoupleTable(count)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
 MSN_TypeI = count(:,1); MSN_TypeII = count(:,2); MSN_TypeIII = count(:,3); MSN_NoIdVTA = count(:,4); 
 FSI_TypeI = count(:,5); FSI_TypeII = count(:,6); FSI_TypeIII = count(:,7); FSI_NoIdVTA = count(:,8);
 CIN_TypeI = count(:,9); CIN_TypeII = count(:,10); CIN_TypeIII = count(:,11); CIN_NoIdVTA = count(:,12);
 NoIdVS_TypeI = count(:,13); NoIdVS_TypeII = count(:,14); NoIdVS_TypeIII = count(:,15); NoIdVS_NoIdVTA = count(:,16);

 TLabels = table(MSN_TypeI, MSN_TypeII, MSN_TypeIII, MSN_NoIdVTA, FSI_TypeI, FSI_TypeII, FSI_TypeIII, FSI_NoIdVTA,...
                 CIN_TypeI, CIN_TypeII, CIN_TypeIII, CIN_NoIdVTA, NoIdVS_TypeI, NoIdVS_TypeII, NoIdVS_TypeIII, NoIdVS_NoIdVTA);

end

