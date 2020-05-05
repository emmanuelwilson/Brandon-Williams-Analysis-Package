%% This Script will convert a graph into a AVI video file
%Simply create your wanted graph and plot it inside the loop on line 12

nframe=100;                                 %# of frames
% x = [ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];       %plot values, this particular instance is for a bar graph, the number of 1s = # of bars
x = [0.1,0.3;0.2,0.5;0.58,0.207;0.4,0.3;0.2,0.8;0.3,0.52;0.5,0.8;0.7,0.2;0.7,0.6;0.25,0.3;0.6,0.48;0.5,0.55;0.4,0.67;0.11,0.61;0.6,0.7]; %series of individual dot coordinates
vid = VideoWriter('DOTs.avi');        %create video file, MODIFY VID NAME HERE
open(vid)                                   %open video file
set(gca,'nextplot','replacechildren')       %Have each plot replace the last
hFig = gcf;                                 %figure/graph data handle
hAx  = gca;                                 %axis handle
for k=1:nframe
    scatter(x(:,1),x(:,2),100000,'filled','k')%plot, scatter plot of large dots
%   bar(x,0.5,'k')                            %plot, INSERT WANTED PLOT HERE 
  axis off                                  %turn off axes 
  % set the figure to full screen
  set(hFig,'units','normalized','outerposition',[0 0 1 1]);
  % set the axes to full screen
  set(hAx,'Unit','normalized','Position',[0 0 1 1]);
  mov=getframe(gcf);                        %save figure into frame
  writeVideo(vid,mov);                      %save frame into video file
end
close(vid)                                  %Stop editing video file