function h = mkLine(d,xl,varargin)    

    % mkLine([{randn(400,10)}; {randn(400,10)}]) <-------- example

    if ~iscell(d)
        d = {d};
%         tmp = repmat({[]},[1 length(d(1,:))]);
%         for i = 1:length(d(:,1))
%             tmp{i} = d(:,i);
%         end
%         d = tmp;
    end
    
     if nargin<2 || isempty(xl)
        xl = 1:length(d{1}(1,:));
    end
    
    groups = length(d);
    ticks = length(d{1}(1,:));
    set(gca,'xlim',[xl(1)-2 xl(end)+2],'xtick',xl,'xticklabel',xl,'fontname','arial',...
        'fontsize',10,'fontweight','bold')
    hold on
    
    w = 0.8;
    wpg = w./groups;
%     groupColor = [{[0.9 0.6 0.6]} {[0.6 0.6 0.9]} {[0.9 0.9 0.6]} {[0.6 0.9 0.6]}];
%     groupColor = [{[0.7 0.7 0.9]} {[0.6 0.6 0.9]} {[0.5 0.5 0.9]} {[0.4 0.4 0.9]} {[0.9 0.9 0.6]} {[0.6 0.9 0.6]} {[0.6 0.9 0.6]} {[0.6 0.9 0.6]}];
%     groupColor = cat(1,groupColor{:});
%     groupColor = (cool(groups)/2)+.3;
%     tmp = transcm;
%     groupColor = tmp(round(linspace(1,length(tmp),groups)),:);
        
    
%     groupColor = [0.5 0.5 0.5; 0.7 0.7 0.7; 0.5 0.5 0.9; 0.7 0.7 0.9; ...
%         0.9 0.5 0.5; 0.9 0.7 0.7; 0.8 0.5 0.8; 0.8 0.7 0.8; 0.5 0.5 0.5; 0.7 0.7 0.7;...
%         [0.5 0.5 0.5; 0.7 0.7 0.7; 0.5 0.5 0.9; 0.7 0.7 0.9; ...
%         0.9 0.5 0.5; 0.9 0.7 0.7; 0.8 0.5 0.8; 0.8 0.7 0.8; 0.5 0.5 0.5; 0.7 0.7 0.7]];
    groupColor = [0.35 0.45 1; 0.85 0.4 1; 0 0 0; 0.75 0.75 0.75; 0.4 0.9 1];
%     groupColor = ones(100,3);
    edgeColor = [0 0 0];
    dotColor = [0.75 0.75 0.75];
%     groupColor = repmat(hsv(4),[4 1]);
%     colors = transcm;
%     groupColor = [colors(([1 40 160 200]),:); 0 0 0 ];

    ebW = 6;
    th = [];
    for i = 1:groups
        doSpot = xl;
        m = nanmean(d{i});
        se = nanstd(d{i})./sqrt(nansum(~isnan(d{i})));
        patch([doSpot(1:length(doSpot)-1); doSpot(1:length(doSpot)-1); ... 
            doSpot(2:length(doSpot)); doSpot(2:length(doSpot))],...
            [m(1:end-1)+se(1:end-1); m(1:end-1)-se(1:end-1); ...
            m(2:end)-se(2:end); m(2:end)+se(2:end)],groupColor(i,:),'edgecolor','none',...
            'facecolor',groupColor(i,:),'edgealpha',0,'facealpha',0.5)
        h(i) = plot(xl,m,'color',groupColor(i,:),'linestyle','-','linewidth',1.25);
%         plot(xl,d{i},'color',0.5+[groupColor(i,:)./2],'linestyle','-','linewidth',1.25)
    end
    
    if nanmin(get(gca,'ylim'))<0
        plot(get(gca,'xlim'),[0 0],'linestyle','--','linewidth',1,'color','k');
    end
end