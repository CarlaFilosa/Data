function [new_data] = parId(new_data,a,nnr)
% ParId re-organize the parameters of the Max's data
% In this file I divided the ID of the units coherently with the
% spiketrains file
% new_ind=newInd;
% nnr=nn;
count=0;
cc1=[];
for i=1:size(a,2)
    
%      if i~=1
%     count=count+nn(a(i-1));
      count=count+nnr(a(i));
      cc1(i)=count;
%     end
%     for j=1:nn(a(i))
     for j=1:nnr(a(i))
         if nnr(a(i))>2
        %nn(a(i))
%         disp('ciao')
        if i==1
            new_data.par{i}(j).ID=new_data.params(j).ID;
            new_data.par{i}(j).animal=new_data.params(j).animal;
            new_data.par{i}(j).date=new_data.params(j).date;
            new_data.par{i}(j).trode=new_data.params(j).trode;
            new_data.par{i}(j).sessionID=new_data.params(j).sessionID;
            new_data.par{i}(j).region=new_data.params(j).region;
            new_data.par{i}(j).ITI_rate=new_data.params(j).ITI_rate;
            new_data.par{i}(j).session_tag=new_data.params(j).session_tag;
            new_data.par{i}(j).LV=new_data.params(j).LV;
            new_data.par{i}(j).CV2=new_data.params(j).CV2;
            new_data.par{i}(j).IR=new_data.params(j).IR;
            new_data.par{i}(j).BP=new_data.params(j).BP;
            new_data.par{i}(j).BT=new_data.params(j).BT;
            new_data.par{i}(j).BE=new_data.params(j).BE;
            new_data.par{i}(j).mean_pause_length=new_data.params(j).mean_pause_length;
            new_data.par{i}(j).mn_pauses_min=new_data.params(j).mn_pauses_min;
            new_data.par{i}(j).fraction_2pauses=new_data.params(j).fraction_2pauses;
            new_data.par{i}(j).bimodality_p=new_data.params(j).bimodality_p;
            new_data.par{i}(j).Fano=new_data.params(j).Fano;
            new_data.par{i}(j).lasertag=new_data.params(j).lasertag;
            new_data.par{i}(j).class=new_data.params(j).class;
            new_data.par{i}(j).stableID=new_data.params(j).stableID;
%             new_data.par{i}(j).safe_region=new_data.params(j).safe_region;
%             new_data.par{i}(j).slowtag=new_data.params(j).slowtag;
%             new_data.par{i}(j).crosstagOT=new_data.params(j).crosstagVTA;
%             new_data.par{i}(j).tetID=new_data.params(j).tetID;
%             new_data.par{i}(j).trode=new_data.params(j).trode;
%             new_data.par{i}(j).sessionID=new_data.params(j).sessionID;
        elseif i~=1
            
            new_data.par{i}(j).ID=new_data.params(j+cc1(i-1)).ID;
            new_data.par{i}(j).animal=new_data.params(j+cc1(i-1)).animal;
            new_data.par{i}(j).date=new_data.params(j+cc1(i-1)).date;
            new_data.par{i}(j).trode=new_data.params(j+cc1(i-1)).trode;
            new_data.par{i}(j).sessionID=new_data.params(j+cc1(i-1)).sessionID;
            new_data.par{i}(j).region=new_data.params(j+cc1(i-1)).region;
            new_data.par{i}(j).ITI_rate=new_data.params(j+cc1(i-1)).ITI_rate;
            new_data.par{i}(j).session_tag=new_data.params(j+cc1(i-1)).session_tag;
            new_data.par{i}(j).LV=new_data.params(j+cc1(i-1)).LV;
            new_data.par{i}(j).CV2=new_data.params(j+cc1(i-1)).CV2;
            new_data.par{i}(j).IR=new_data.params(j+cc1(i-1)).IR;
            new_data.par{i}(j).BP=new_data.params(j+cc1(i-1)).BP;
            new_data.par{i}(j).BT=new_data.params(j+cc1(i-1)).BT;
            new_data.par{i}(j).BE=new_data.params(j+cc1(i-1)).BE;
            new_data.par{i}(j).mean_pause_length=new_data.params(j+cc1(i-1)).mean_pause_length;
            new_data.par{i}(j).mn_pauses_min=new_data.params(j+cc1(i-1)).mn_pauses_min;
            new_data.par{i}(j).fraction_2pauses=new_data.params(j+cc1(i-1)).fraction_2pauses;
            new_data.par{i}(j).bimodality_p=new_data.params(j+cc1(i-1)).bimodality_p;
            new_data.par{i}(j).Fano=new_data.params(j+cc1(i-1)).Fano;
            new_data.par{i}(j).lasertag=new_data.params(j+cc1(i-1)).lasertag;
            new_data.par{i}(j).class=new_data.params(j+cc1(i-1)).class;
            new_data.par{i}(j).stableID=new_data.params(j+cc1(i-1)).stableID;
%             new_data.par{i}(j).safe_region=new_data.params(j+cc1(i-1)).safe_region;
%             new_data.par{i}(j).slowtag=new_data.params(j+cc1(i-1)).slowtag;
%             new_data.par{i}(j).crosstagOT=new_data.params(j+cc1(i-1)).crosstagOT;
%             new_data.par{i}(j).crosstagVTA=new_data.params(j+cc1(i-1)).crosstagVTA;
%             new_data.par{i}(j).tetID=new_data.params(j+cc1(i-1)).tetID;
        end
         end
    end
    
end
end

