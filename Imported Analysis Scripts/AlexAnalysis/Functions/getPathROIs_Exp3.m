function getPathROIs(paths)
    clc
    fprintf(['Define Path ROIs:\n']);
    
    doorROI = nan(4,3);
    roomROI = nan(4,1);
    for p = paths'
        s = load(p{1});
        
%         if isfield(s.processed,'roi')
%             continue
%         end
        
        fprintf(['\t' num2str(p{1}) '\n'])
        
        figure(1)
        set(gcf,'position',[50 50 nanmax(s.processed.p').*10])
        plot(s.processed.p(1,:),s.processed.p(2,:));
        hold on
        set(gca,'xlim',[0 nanmax(s.processed.p(1,:))],'ylim',[0 nanmax(s.processed.p(2,:))])
        drawnow
        
        fprintf(['\t\tSelect Room.'])
        rect = getrect(); 
        roomROI(:,1) = rect; %uint16([rect(1) rect(1)+rect(3) rect(2) rect(2)+rect(4)]);
        
        fprintf(['\tSelect Left Door.'])
        rect = getrect(); 
        doorROI(:,1) = rect; %uint16([rect(1) rect(1)+rect(3) rect(2) rect(2)+rect(4)]);
        
        fprintf(['\tSelect Top Door.'])
        rect = getrect(); 
        doorROI(:,2) = rect; %uint16([rect(1) rect(1)+rect(3) rect(2) rect(2)+rect(4)]);
        
        fprintf(['\tSelect Long Entrance Door.'])
        rect = getrect(); 
        doorROI(:,3) = rect; %uint16([rect(1) rect(1)+rect(3) rect(2) rect(2)+rect(4)]);
        
%         fprintf(['\tSelect Exclusion Breach.'])
%         rect = getrect(); 
%         doorROI(:,4) = rect; %uint16([rect(1) rect(1)+rect(3) rect(2) rect(2)+rect(4)]);
        
%         save('PathROIs','doorROI','roomROI');
%         
%         load('PathROIs','doorROI','roomROI');
        
        s.processed.roi.room = roomROI;
        s.processed.roi.door = doorROI;
        
        for i = 1:length(roomROI(1,:))
            rectangle('position',roomROI(:,i),'linewidth',2);
        end
        for i = 1:length(doorROI(1,:))
            rectangle('position',doorROI(:,i),'linewidth',2,'linestyle','--');
        end
        drawnow
        pause(5)
        close all
        drawnow
        
        save(p{1},'-struct','s','-v7.3');
    end
end