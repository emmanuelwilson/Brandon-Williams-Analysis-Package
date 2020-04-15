function cleanPath(paths)
    clc
    fprintf(['Clean Path:\n']);
    
    doorROI = nan(4,6);
    roomROI = nan(4,3);
    for p = paths'
        s = load(p{1});
        
        fprintf(['\t' num2str(p{1}) '\n'])
        
        userInput = 'Y';  
        
        unP = s.pos.uninterp([1 2],:);
        unP = bsxfun(@minus,unP,nanmin(unP')');
        unP = bsxfun(@times,unP,nanmax(s.environment.size./nanmax(unP')));
        oldP = unP;
        
        linX = linterp(s.processed.posFrames(:,2),oldP(1,s.processed.posFrames(:,1)),s.processed.validTraceFrames(:,2));
        linY = linterp(s.processed.posFrames(:,2),oldP(2,s.processed.posFrames(:,1)),s.processed.validTraceFrames(:,2));

        unX = linterp(s.processed.posFrames(:,2),unP(1,s.processed.posFrames(:,1)),s.processed.validTraceFrames(:,2));
        unY = linterp(s.processed.posFrames(:,2),unP(2,s.processed.posFrames(:,1)),s.processed.validTraceFrames(:,2));
        s.processed.p = [unX'; unY'];

        %%%% For Ephys Data
        
%         unP = s.pos.p([1 2],:);
%         unP = bsxfun(@minus,unP,nanmin(unP')');
%         unP = bsxfun(@times,unP,((s.environment.size)./nanmax(unP'))');
%         s.processed.p = unP;
        
        %%%%%%

        
%         s.processed.p = bsxfun(@minus,s.processed.p,[s.environment.size./2]');
%         while (strcmp(userInput,'Y')) 
%             hold off
%             figure(1)
%             set(gcf,'position',[50 50 nanmax(s.processed.p').*20])
%             plot(s.processed.p(1,:),s.processed.p(2,:),'linestyle','none','marker','o',...
%                 'markersize',2,'color','k');
%             hold off
%             set(gca,'xlim',[nanmin(s.processed.p(1,:)) nanmax(s.processed.p(1,:))],...
%                 'ylim',[nanmin(s.processed.p(2,:)) nanmax(s.processed.p(2,:))])
%             drawnow
%             
%             str = input('rotate?','s');
%             
%             if  str == 'r'
%                 s.processed.p = [cosd(1) -sind(1); sind(1) cosd(1)]*s.processed.p;
%             elseif str == 'l'
%                 s.processed.p = [cosd(-1) -sind(-1); sind(-1) cosd(-1)]*s.processed.p;
%             else
%                 userInput = 'N';
%                 continue
%             end
%             
%             
%         end
        
        userInput = 'Y';  
        [s.processed.p gaps] = interpNaNs(s.processed.p');
        s.processed.p = s.processed.p';
%         gapSize = 60;
%         while any(gaps>gapSize)
%             start = find(gaps>gapSize,1,'first');
%             stop = find(gaps(start:end)<gapSize,1,'first')-1;
%             if isempty(stop)
%                 stop = length(gaps)-start;
%             end
%             
%             hold off
%             figure(1)
%             set(gcf,'position',[50 50 nanmax(s.processed.p').*20])
%             plot(s.processed.p(1,:),s.processed.p(2,:),'linestyle','none','marker','o',...
%                 'markersize',2,'color',[0.6 0.6 0.6]);
%             hold on
%             plot(s.processed.p(1,start:start+stop),s.processed.p(2,start:start+stop),'linestyle','none','marker','o',...
%                 'markersize',3,'color','k');
%             hold off
%             set(gca,'xlim',[nanmin(s.processed.p(1,:)) nanmax(s.processed.p(1,:))],...
%                 'ylim',[nanmin(s.processed.p(2,:)) nanmax(s.processed.p(2,:))])
%             drawnow
%             
%             [ax ay] = getpts;    
%             
%             if length(ax)>=2
%                 
%                 inds = knnsearch(s.processed.p(:,start:start+stop)',[ax ay]);
%                 
%                 s.processed.p(:,start:start+stop-1) = nan;
%                 for i = 1:2:length(inds)
%                     s.processed.p(:,start+inds(i)-1) = [ax(i+1) ay(i+1)];
%                 end
%                 
%                 s.processed.p = interpNaNs(s.processed.p')';
%             end  
%             gaps(start:start+stop) = 0;
%         end
        
        while (strcmp(userInput,'Y')) 
            hold off
            figure(1)
            set(gcf,'position',[50 50 nanmax(s.processed.p').*20])
            plot(s.processed.p(1,:),s.processed.p(2,:),'linestyle','none','marker','o',...
                'markersize',2,'color','k');
            hold off
            set(gca,'xlim',[nanmin(s.processed.p(1,:)) nanmax(s.processed.p(1,:))],...
                'ylim',[nanmin(s.processed.p(2,:)) nanmax(s.processed.p(2,:))])
            drawnow
            
            [ax ay] = getpts;    
            
            if length(ax)~=4
                userInput = 'N';
                continue
            end
            
            inds = knnsearch(s.processed.p',[ax ay]);
            
            s.processed.p(:,inds(3)) = s.processed.p(:,inds(4));
            s.processed.p(:,[min(inds(1),inds(2)):inds(3)-1 ...
                inds(3)+1:max(inds(1),inds(2))]) = nan;
            s.processed.p = interpNaNs(s.processed.p')';
            
        end
        
        userInput = 'Y';  
        
        while (strcmp(userInput,'Y')) 
            figure(1)
            plot(s.processed.p(1,:),s.processed.p(2,:),'linestyle','none','marker','o',...
                'markersize',2,'color','k');
            hold off
            set(gca,'xlim',[nanmin(s.processed.p(1,:)) nanmax(s.processed.p(1,:))],...
                'ylim',[nanmin(s.processed.p(2,:)) nanmax(s.processed.p(2,:))])
            drawnow
            
            rect = getrect();
            exclude = s.processed.p(1,:)>rect(1) & s.processed.p(1,:) < rect(1)+rect(3) & ...
                s.processed.p(2,:)>rect(2) & s.processed.p(2,:) < rect(2)+rect(4);
            
            
            s.processed.p(:,exclude) = nan;
            [s.processed.p gaps] = interpNaNs(s.processed.p');
            s.processed.p = s.processed.p';
%             s.processed.p(:,gaps>15) = nan;
            
            if ~any(exclude)
                userInput = 'N';
            end
        end
        
        close all
        drawnow
        
        s.processed.p = bsxfun(@minus,s.processed.p,nanmin(s.processed.p')');
        s.processed.p = bsxfun(@rdivide,s.processed.p,nanmax(s.processed.p')');
        s.processed.p = bsxfun(@times,s.processed.p,s.environment.size');
        
        save(p{1},'-struct','s','-v7.3');
    end
end