function h = mkGraph(d,xl,doDots,varargin)
    h = [];
    if nargin <3 | isempty(doDots)
        doDots = true;
    end
    
    if nargin<2 || isempty(xl)
        xl = 1:length(d(1,:));
    end
    
    if ~iscell(d)
        tmp = repmat({[]},[1 length(d(1,:))]);
        for i = 1:length(d(1,:))
            tmp{i} = d(:,i);
        end
        d = tmp;
    end
    
    groups = length(d(:,1));
    ticks = length(d(1,:));
    set(gca,'xlim',[0.5 ticks+0.5],'xtick',[1:ticks],'xticklabel',xl,'fontname','arial',...
        'fontsize',10,'fontweight','bold')
    hold on
    
    w = 0.4;
    wpg = w./groups;
%     groupColor = [{[0.9 0.6 0.6]} {[0.6 0.6 0.9]} {[0.9 0.9 0.6]} {[0.6 0.9 0.6]}];
%     groupColor = [{[0.7 0.7 0.9]} {[0.6 0.6 0.9]} {[0.5 0.5 0.9]} {[0.4 0.4 0.9]} {[0.9 0.9 0.6]} {[0.6 0.9 0.6]} {[0.6 0.9 0.6]} {[0.6 0.9 0.6]}];
%     groupColor = cat(1,groupColor{:});
%     groupColor = (cool(groups)/2)+.3;
%     tmp = transcm;
%     groupColor = tmp(round(linspace(1,length(tmp),groups)),:);
        
    
    groupColor = [0.9 0.3 0.3; 0.3 0.3 0.9; 0.2 0.2 0.2; 0.5 0.5 0.9; 0.7 0.7 0.9; ...
        0.9 0.5 0.5; 0.9 0.7 0.7; 0.8 0.5 0.8; 0.8 0.7 0.8; 0.5 0.5 0.5; 0.7 0.7 0.7;...
        [0.5 0.5 0.5; 0.7 0.7 0.7; 0.5 0.5 0.9; 0.7 0.7 0.9; ...
        0.9 0.5 0.5; 0.9 0.7 0.7; 0.8 0.5 0.8; 0.8 0.7 0.8; 0.5 0.5 0.5; 0.7 0.7 0.7]];
%     groupColor = [0.5 0.5 0.5; 0.5 0.5 0.5; 0.0 0.9 1; 0.75 0.75 0.75; 0.4 0.9 1];
%     groupColor = ones(100,3);
    edgeColor = [1 0.2 0.2; 0.2 0.2 1; 0.2 1 0.2; 0.2 1 1];
    dotColor = [0.9 0.7 0.7; 0.7 0.7 0.9; 0.7 0.9 0.7; 0.7 0.9 0.9];
%     groupColor = repmat(hsv(4),[4 1]);
%     colors = transcm;
%     groupColor = [colors(([1 40 160 200]),:); 0 0 0 ];

    if length(d(:,1))==1
        groupColor = [0 0 0];
        dotColor = [0.2 0.2 0.2];
    end

    ebW = 6;
    th = [];
    jitters = repmat({[]},[size(d)]);
    for i = 1:groups
        for j = 1:ticks
            jitters{i,j} = randn(size(d{i,j})).*0.035;
        end
    end
    for i = 1:groups
        for j = 1:ticks
            
            if length(d{i,j})>=1
                if doDots
                    plot(j-(w./2)+(wpg.*(i-1))+wpg./2+jitters{i,j},d{i,j},'linestyle','none','color',dotColor(i,:),...
                        'marker','o','markersize',5);
                end
                if j < ticks && length(d{i,j}) == length(d{i,j+1}) && length(d(:,1)) == 1
%                     plot([repmat([j j+1],[length(d{i,j}) 1])+[jitters{i,j} jitters{i,j+1}]]',[d{i,j} d{i,j+1}]', ...
%                         'color',dotColor(i,:),'linewidth',0.8)
                end
                if i < groups && length(d{i,j}) == length(d{i+1,j}) && length(d(:,1)) > 1
%                     plot([repmat([j-(w./2)+(wpg.*(i-1))+wpg./2 ...
%                         j-(w./2)+(wpg.*(i))+wpg./2],[length(d{i,j}(:,1)) 1])+[jitters{i,j} jitters{i+1,j}]]', ...
%                         [d{i,j} d{i+1,j}]','color',dotColor(i,:),'linewidth',0.8)
                end
            end
            
        end
    end
    
    for i = 1:groups
        gm = cellfun(@nanmedian,d(i,:));
        j = [1:ticks]';
%         xs = nanmean([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i))],2)';
        h(i) = plot(j,gm,'color',edgeColor(i,:)','linewidth',1);
        for j = 1:ticks
            
            dir = 1;
%             h(i,j) = patch([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i)) (j-(w./2))+(wpg.*(i)) (j-(w./2))+(wpg.*(i-1))],...
%                 [0 0 repmat(nanmean(d{i,j}),[1 2])],groupColor(ceil(j),:)./(i),'edgecolor',edgeColor,'linewidth',1.5);
%             set(h(i,j),'facealpha',0.5);
% % %             plot([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i))],repmat(nanmean(d{i,j}),[1 2]),...
% % %                 'color',edgeColor(i,:)','linewidth',1)
% % %             if length(d{i,j})>1
% % %                 plot(repmat(nanmean([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i))]),[1 2]),...
% % %                     [nanmean(d{i,j})-dir.*nanstd(d{i,j})./sqrt(sum(~isnan(d{i,j}))) ...
% % %                     nanmean(d{i,j})+dir.*nanstd(d{i,j})./sqrt(sum(~isnan(d{i,j})))], ...
% % %                     'color',edgeColor(i,:),'linewidth',1);
% % %             end
            
            if ~isempty(varargin) && ~isempty(varargin{1})
                th(end+1) = text(nanmean([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i)) ...
                    (j-(w./2))+(wpg.*(i)) (j-(w./2))+(wpg.*(i-1))]),...
                    0,[' ' varargin{1}{i}],'horizontalalignment','left',...
                    'verticalalignment','middle','fontname','arial','fontsize',9,'rotation',90,...
                    'fontweight','bold','color','k');
            end
        end
    end
    
%     
%     scale = range(get(gca,'ylim')).*0.05;
%     maxd = cellfun(@nanmean,d)+cellfun(@nanstd,d)./sqrt(cellfun(@numel,d))+scale;
%     numBars = ones(size(d));
%     for i = 1:groups
%         for j = 1:ticks
%             for ii = groups:-1:i
%                 for jj = j:ticks
%                     try [hyp p] = ttest2(d{i,j},d{ii,jj},0.05,'both');
%                         if p<0.05
%                             ind1 = sub2ind(size(d),i,j);
%                             ind2 = sub2ind(size(d),ii,jj);
%                             currMax = max(maxd(ind1:ind2));
%                             h2 = plot([nanmean([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i))]) ...
%                                 nanmean([(jj-(w./2))+(wpg.*(ii-1)) (jj-(w./2))+(wpg.*(ii))])],...
%                                 ones(1,2).*currMax,'color','k');
%                             plot([repmat(nanmean([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i))]),[1 2])],...
%                                 [currMax currMax-scale./2],'color','k')
%                             plot([repmat(nanmean([(jj-(w./2))+(wpg.*(ii-1)) (jj-(w./2))+(wpg.*(ii))]),[1 2])],...
%                                 [currMax currMax-scale./2],'color','k')
%                             maxd(ind1:ind2) = max(maxd(ind1:ind2))+scale;
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     dh = get(gca,'Children');
%     ys = [];
%     for i = 1:length(dh)
%         try ys = [ys; get(dh(i),'YData')];
%         catch
%         end
%     end
%     set(gca,'ylim',[min(get(gca,'ylim')) ...
%         max(ys(:))+scale]);
%     if min(get(gca,'ylim'))<0
%         plot(get(gca,'xlim'),[0 0],'color','k','linewidth',1)
%     end
    scale = range(get(gca,'ylim')).*0.05;
    for i = 1:groups
        for ii = 1+1:groups
            for j = 1:ticks
                try [hyp p] = ttest2(d{i,j},d{ii,j},0.05,'both');
%                     if p<0.05
%                         ind = max([nanmean(d{i,j})+((mod(i,2).*2)-1).*nanstd(d{i,j})./(sqrt(length(d{i,j}))),...
%                             nanmean(d{ii,j})+((mod(ii,2).*2)-1).*nanstd(d{ii,j})./(sqrt(length(d{ii,j}))),...
%                             nanmean(d{i,j}),nanmean(d{ii,j})]).*1+scale;
%                         text(j,ind,'*','fontname','arial','fontsize',12,'fontweight','bold','horizontalalignment','center',...
%                             'verticalalignment','middle')
%                     end
                end
            end
        end
    end
    
    if nanmin(get(gca,'ylim'))<0
        plot(get(gca,'xlim'),[0 0],'linestyle','--','linewidth',1,'color','k');
    end
end