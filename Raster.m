function [myplot] = Raster(SpM,color_par,STR)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% figure(1)
% if nargin <2 || isempty(color_par) color_par=[0 0 1]; end

for  j= 1:size(SpM,1)
    box on;
myplot=plot(SpM(j,:),j*ones(size(SpM(j,:))),'b.')
% plot(spiketrain{i}(j,:),repmat(j,size(spiketrain{i}(j,:))),'b.')
% myplot=plot(SpM(j,:),repmat(j,size(SpM(j,:))),'.','color',color_par)
if nargin==3
    title(STR)
end
hold on
xlim([min(min(SpM)),max(max(SpM))])
ylim([1,size(SpM,1)])
end

end

