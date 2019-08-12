load('c_Stable')
VTA_Laser=c_Stable.VTA_ltag;
VTA_LaserIn=c_Stable.VTA_linh;
type2=c_Stable.type2;
gg=0;
for i=1:length(VTA_Laser)
    for j=1:length(type2)
       if VTA_Laser(i)==type2(j)
           gg=gg+1;
           c_Stable.Type2L(gg)=type2(j);
       end
    end
end

gg=0;
for i=1:length(VTA_LaserIn)
    for j=1:length(type2)
       if VTA_LaserIn(i)==type2(j)
           gg=gg+1;
           c_Stable.Type2LIn(gg)=type2(j);
       end
    end
end

clear gg i j