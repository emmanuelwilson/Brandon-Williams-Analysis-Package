
% try
%     load('SocialProximity_V2.mat')
% catch
%     try
        load('CupStats.mat')
        load('ms.mat');
        load('HeadTrackingData.mat');
        load('frameMap.mat');
        load('badframes.mat')
        frameMaptemp = frameMap(t:end);       
        
% %         h = figure(1);
% %         vid = VideoReader('behavCam1.avi');
% %         f1 = vid.read(1);
% %         h = figure;
% %         imshow(f1);
% %         hold on
% %         plot(SINKdata(frameMaptemp,1),SINKdata(frameMaptemp,2))
% %         scatter(out.QPOL(1),out.QPOL(2))
% %         scatter(out.QPOR(1),out.QPOR(2))
% %         scatter(out.QPW(1,1),out.QPW(1,2))
% %         scatter(out.QPW(2,1),out.QPW(2,2))
% %         scatter(out.QPW(3,1),out.QPW(3,2))
% %         scatter(out.QPW(4,1),out.QPW(4,2))
% %         [ObjectPos1, ObjectPos2] = getpts(h);
% %         Objects = cat(2,ObjectPos1,ObjectPos2);
% %         Object1 = Objects(1,:);
% %         Object2 = Objects(2,:);
% %         [cornerPos1, cornerPos2] = getpts(h);
% %         QPW = cat(2,cornerPos1, cornerPos2);
% %         close(h);
% %         
% %         out.QPW = QPW;
% %         out.QPOL= Object1;
% %         out.QPOR= Object2;
% %         save('CupStats.mat','out')
% %         
%     catch       
%         load('HeadTrackingData.mat');
%         load('frameMap.mat');
%         load('badframes.mat')        
%         frameMaptemp = frameMap(t:end);
%         QPW = findEdges(SINKdata(frameMaptemp,:));
%         h = figure;
%         plot(SINKdata(frameMaptemp,1),SINKdata(frameMaptemp,2))
%         [ObjectPos1, ObjectPos2] = getpts(h);
%         Objects = cat(2,ObjectPos1,ObjectPos2);
%         Object1 = Objects(1,:);
%         Object2 = Objects(2,:);
%         close(h);
%         
%         out.QPW = QPW;
%         out.QPOL= Object1;
%         out.QPOR= Object2;
%         save('CupStats.mat','out')
%     end   
   SocialProximity = SocialProximityFiring_V3(ms,SINKdata,frameMap,out, 5, 3.75, 45, t);
   save('SocialProximity_V2.mat','SocialProximity')
    
% end


function QP = findEdges(tracking)
ifEscape = 0;
h=figure(1);

while ~ifEscape
    figure(h);
    clf
    
    %[occupancy, xdim, ydim]=root.Occupancy([],[],1,2);
    %imagesc(xdim,ydim,occupancy');
    set(gca,'YDir','Normal'); %colormap(jet);
    clim=get(gca,'clim');set(gca,'clim',clim/50);
    hold on
    plot(tracking(:,1),tracking(:,2),'k');
    QP = [];
    
    set(h,'Name','Select Corners of Walls. Esc--> done. **Do not complete!**')
    
    button = 1;
    
    while button~=27
        [x,y,button] = ginput(1);
        
        clf
        
        %imagesc(xdim,ydim,occupancy');
        set(gca,'YDir','Normal'); %colormap(jet);
        clim=get(gca,'clim');set(gca,'clim',clim/50);
        hold on
        plot(tracking(:,1),tracking(:,2),'k');
        
        if ~isempty(QP)
            plot(QP(:,1),QP(:,2),'r')
            plot(QP(:,1),QP(:,2),'ro','MarkerFaceColor','r')
        end
        
        if button == 32 %space bar
            QP = [QP; NaN NaN];
        elseif button~=27
            QP = [QP; x y];
        end
        
        plot(QP(:,1),QP(:,2),'r')
        plot(QP(:,1),QP(:,2),'ro','MarkerFaceColor','r')
        
    end
    
    %Ask for verification
    edg = splitter(QP);
    clf;
    set(h,'Name','Verify. 0--> Try again; 1--> Confirm')
    plot(tracking(:,1),tracking(:,2),'k');
    hold on
    
    for m = 1:numel(edg)
        for n = 1:size(edg{m},1)
            sp = squeeze(edg{m}(n,:,1));
            ep = squeeze(edg{m}(n,:,2));
            plot([sp(1) ep(1)],[sp(2) ep(2)],'ro','MarkerFaceColor','r')
            plot([sp(1) ep(1)],[sp(2) ep(2)],'r')
        end
    end   
    % set or repeat
    while button ~=48 && button~=49
        [~,~,button]=ginput(1);
    end
    ifEscape = button==49;
    
end

close(h);
drawnow();
end

%Split the corner coordinates in X and Y vectors
function edg = splitter(QP)

inds = find(isnan(QP(:,1)));
xs=CMBHOME.Utils.SplitVec(QP(:,1), @(x) isnan(x));
ys=CMBHOME.Utils.SplitVec(QP(:,2), @(x) isnan(x));

% split corners
for m = 1:size(xs,1)
    QP2{m} = [xs{m} ys{m}];
    QP2{m}(find(isnan(QP2{m}(:,1))),:) = [];
end

for m = 1:numel(QP2)
    for n = 1:size(QP2{m},1)
        sp = n;ep=n+1;
        if ep>size(QP2{m},1), ep=1;end
        edg{m}(n,:,1) = [QP2{m}(sp,1) QP2{m}(sp,2)];
        edg{m}(n,:,2) = [QP2{m}(ep,1) QP2{m}(ep,2)];
    end
end

end