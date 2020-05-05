function aggPaths(paths)
    piece = [];
    spiece = [];
    
    
    warning off all
    
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'))-1;
        piece = [piece; {paths{i}(ind(end-1)+2:ind(end))}];
    end
    upiece = unique(piece);
    
    order = [{'A'} {'A2'} {'B'} {'A3'}];
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}) '\n'])
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);
        figure
        set(gcf,'position',[50 50 250.*(length(sessions)) 250])
        
        for si = 1:length(sessions)
            s = load(sessions{si},'processed','properties');
            subplot(2,length(sessions),find(ismember(order,s.properties.session)))
            plot(s.processed.p(1,:),s.processed.p(2,:),'color','k','linewidth',0.75,...
                'linestyle','-')
            hold on
%             doC = hsv(length(s.env.nlocs(1,:)));
%             for j = 1:length(s.env.nlocs(1,:))
%                 hi(j) = plot(s.env.nlocs(1,j),s.env.nlocs(2,j),'marker','o','linestyle','none',...
%                     'markersize',4,'markeredgecolor',doC(j,:),'markerfacecolor',doC(j,:));
%             end
            rectangle('position',[0 0 s.processed.envSize],'edgecolor',[0.5 0.5 0.5])
            set(gca,'xlim',[-1 75+1],'ylim',[-1 75+1])
%             legend(hi,cellfun(@num2str,num2cell(1:length(s.env.nlocs(1,:))),'uniformoutput',false), ...
%                 'location','southoutside','orientation','horizontal')
            axis square
            axis off
        end
        outP = ['Plots/Agg/Trajectories/' slind(sessions{1},[2 0])];
        saveFig(gcf,outP,[{'tiff'} {'pdf'}])
        close all
        drawnow
    end
end