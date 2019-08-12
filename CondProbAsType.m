function [TabCond,CondProb] = CondProbAsType(JointProb, AsNeuTypeProb,given)
% I calculate here the conditional probability P(AS_typeN_{k}|AS_typeN_{i})


CondProb = JointProb./AsNeuTypeProb;

switch given
 case 'vs|vta'
MSN_g_TypeI = CondProb(:,1); MSN_g_TypeII = CondProb(:,2); MSN_g_TypeIII = CondProb(:,3); MSN_g_NoIdVTA= CondProb(:,4);
FSI_g_TypeI = CondProb(:,5); FSI_g_TypeII = CondProb(:,6); FSI_g_TypeIII = CondProb(:,7); FSI_g_NoIdVTA= CondProb(:,8);
CIN_g_TypeI = CondProb(:,9); CIN_g_TypeII = CondProb(:,10); CIN_g_TypeIII = CondProb(:,11); CIN_g_NoIdVTA= CondProb(:,12);
NoIdVS_g_TypeI = CondProb(:,13); NoIdVS_g_TypeII = CondProb(:,14); NoIdVS_g_TypeIII = CondProb(:,15); NoIdVS_g_NoIdVTA= CondProb(:,16);

TabCond = table(MSN_g_TypeI,MSN_g_TypeII,MSN_g_TypeIII,MSN_g_NoIdVTA,FSI_g_TypeI,FSI_g_TypeII,FSI_g_TypeIII,FSI_g_NoIdVTA,...
CIN_g_TypeI,CIN_g_TypeII,CIN_g_TypeIII,CIN_g_NoIdVTA, NoIdVS_g_TypeI,NoIdVS_g_TypeII,NoIdVS_g_TypeIII,NoIdVS_g_NoIdVTA);

case 'vta|vs'

TypeI_g_MSN = CondProb(:,1); TypeI_g_FSI = CondProb(:,5); TypeI_g_CIN = CondProb(:,9); TypeI_g_NoIdVS = CondProb(:,13);
TypeII_g_MSN = CondProb(:,2); TypeII_g_FSI = CondProb(:,6); TypeII_g_CIN = CondProb(:,10); TypeII_g_NoIdVS = CondProb(:,14);
TypeIII_g_MSN = CondProb(:,3); TypeIII_g_FSI = CondProb(:,7); TypeIII_g_CIN = CondProb(:,11); TypeIII_g_NoIdVS = CondProb(:,15);
NoIdVTA_g_MSN = CondProb(:,4); NoIdVTA_g_FSI = CondProb(:,8); NoIdVTA_g_CIN = CondProb(:,12); NoIdVTA_g_NoIdVS = CondProb(:,16);

TabCond = table(TypeI_g_MSN, TypeII_g_MSN, TypeIII_g_MSN, NoIdVTA_g_MSN, TypeI_g_FSI, TypeII_g_FSI, TypeIII_g_FSI, NoIdVTA_g_FSI,...
TypeI_g_CIN,TypeII_g_CIN,TypeIII_g_CIN,NoIdVTA_g_CIN,TypeI_g_NoIdVS,TypeII_g_NoIdVS,TypeIII_g_NoIdVS,NoIdVTA_g_NoIdVS);

end

end

