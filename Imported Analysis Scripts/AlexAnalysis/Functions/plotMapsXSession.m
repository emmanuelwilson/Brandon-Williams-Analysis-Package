function plotMapsXSession(folder)
    mice = dir(folder);
    mice = {mice(3:end).name};
    
    similarity = repmat({[]},[1 length(mice)]);
    for m = mice
        mp = [folder '/' m{1}];
        sessions = dir(mp);
        sessions = {sessions(3:end).name};
        
        % order sessions by date
        clear dates
        for j = 1:length(sessions)
            dates(j) = datetime(str2num(sessions{j}(1:2)),str2num(sessions{j}(4:5)),str2num(sessions{j}(7:8)));
        end
        [a b] = sort(dates);
        sessions = sessions(b);
        
        batchSize = 6; 
        
        allM = repmat({[]},[1 length(sessions)]);
        allGT = repmat({[]},[1 length(sessions)]);
        allP = repmat({[]},[1 length(sessions)]);
        for dci = 1:length(sessions)
            s = load([mp '/' sessions{dci}]);
            allM{dci} = mkTraceMaps(s.processed.p,s.processed.trace(s.processed.isAligned,:),[]);
            allGT{dci} = s.processed.trace(s.processed.isAligned,:);
            allP{dci} = s.processed.p;
        end
        
        for k = 1:length(allM{1}(1,1,:))
            figure(1)
            set(gcf,'position',[50 50 400.*length(sessions) 200])
            for dci = 1:length(sessions)
                subplot(1,length(sessions).*2,(dci-1).*2 + 1)
                imagesc(allM{dci}(:,:,k))
                set(gca,'ydir','normal')
                colormap('jet')
                alpha(double(~isnan(allM{dci}(:,:,k))))
                axis tight; axis equal; axis off;
                
                
                subplot(1,length(sessions).*2,(dci-1).*2 + 2)
                plot(allP{dci}(2,:),allP{dci}(1,:),'color',[0.6 0.6 0.6],'linewidth',1)
                hold on
                plot(allP{dci}(2,logical(allGT{dci}(k,:))), ...
                    allP{dci}(1,logical(allGT{dci}(k,:))),'color',[0.9 0.2 0.2],...
                    'linestyle','none','marker','o','markerfacecolor',[0.9 0.2 0.2],'markersize',1);
                axis tight; axis equal; axis off;
                drawnow
            end
            
            outP = ['Plots/MapsXSessions/' m{1} '/Cell_' num2str(k)];
            saveFig(gcf,outP,'tiff')
            close all
            drawnow
        end
        
    end
end