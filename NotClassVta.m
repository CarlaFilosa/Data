% Not_Class_VTA = [105;107;154;184;185;186;194;212;796;797;798;800;831;834;854;898;899;948;953;962;963;964;984;985;1029;1379;1403;1417;1418;1419;1444;1445;1468;1530;1551;1552;1563;1583;1584;1586;1629;1633]
% sort(Not_Class_VTA)
% Not_Class_VTA=ans;

save Not_Class_VTA.mat Not_Class_VTA
Not_Class_Ind = zeros(1,length(Not_Class_Vta))
length(d.clust_params)
ll=0;
for i =1:length(d.clust_params)
    for j=1:length(Not_Class_Vta)
    if d.clust_params(i).stableID==Not_Class_Vta(j)
        ll=ll+1;
        Not_Class_Ind(ll)=i;
    end
    end
end


Not_Class_Ind=Not_Class_Ind'
for i=1:length(Not_Class_Ind)
regionNonClass(i) = d.clust_params(Not_Class_Ind).region
end

regionNonClass=regionNonClass'

save NotClass_VTA_Index_NONstable.mat  Not_Class_Ind regionNonClass


rt1=[];
rt2=[];
rt3=[];
oo=0;
ll=0;
kk=0;
for i=1:length(c_Stable.type1)
    for j=1:length(Not_Class_VTA)
        if c_Stable.type1(i)==Not_Class_VTA(j)
             oo=oo+1;
           rt1(oo)=c_Stable.type1(i);
        end
    end
end
for i=1:length(c_Stable.type2)
    for j=1:length(Not_Class_VTA)
     if c_Stable.type2(i)==Not_Class_VTA(j)
             ll=ll+1;
           rt2(ll)=c_Stable.type2(i);
     end
    end
end
           
for i=1:length(c_Stable.type3)
    for j=1:length(Not_Class_VTA)
           if c_Stable.type3(i)==Not_Class_VTA(j)
             kk=kk+1;
           rt3(kk)=c_Stable.type3(i);
           end
    end
end