function h = mkWhisker(d,xl,varargin)
    if nargin<2 || isempty(xl)
        xl = 1:length(d(1,:));
    end

    if ~iscell(d)
        newD = repmat({[]},[1 length(d(1,:))]);
        for i = 1:length(d(1,:))
            newD{i} = d(:,i);
        end
        d = newD;
    end
    
    groups = length(d(:,1));
    ticks = length(d(1,:));
    set(gca,'xlim',[0.5 ticks+0.5],'xtick',[1:ticks],'xticklabel',xl,'fontname','arial',...
        'fontsize',10,'fontweight','bold')
    hold on
    
    w = 0.5;
    wpg = w./groups;
%     groupColor = ones(groups,3).*0.8
    groupColor = [0.9 0.5 0.5; 0.5 0.5 0.9];
%     groupColor = [{[0.7 0.7 0.9]} {[0.6 0.6 0.9]} {[0.5 0.5 0.9]} {[0.4 0.4 0.9]} {[0.9 0.9 0.6]} {[0.6 0.9 0.6]} {[0.6 0.9 0.6]} {[0.6 0.9 0.6]}];
%     groupColor = cat(1,groupColor{:});
%     groupColor = (cool(groups)/2)+.3;
%     tmp = transcm;
%     groupColor = tmp(round(linspace(1,length(tmp),groups)),:);
    
%     groupColor = [0.5 0.5 0.5; 0.7 0.7 0.7; 0.5 0.5 0.9; 0.7 0.7 0.9; ...
%         0.9 0.5 0.5; 0.9 0.7 0.7; 0.8 0.5 0.8; 0.8 0.7 0.8; 0.5 0.5 0.5; 0.7 0.7 0.7;...
%         [0.5 0.5 0.5; 0.7 0.7 0.7; 0.5 0.5 0.9; 0.7 0.7 0.9; ...
%         0.9 0.5 0.5; 0.9 0.7 0.7; 0.8 0.5 0.8; 0.8 0.7 0.8; 0.5 0.5 0.5; 0.7 0.7 0.7]];
    
%     groupColor = repmat(hsv(4),[4 1]);

    ebW = 6;
    th = [];
    for i = 1:groups
        for j = 1:ticks
            if isempty(d{i,j})
                continue
            end
            sVals = sort(d{i,j});
            sVals(isnan(sVals)) = [];
            
%             plot(mean([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i))]).*ones(1,2),...
%                 sVals([1 end]),'color','k','linewidth',2);

            % 1 SD Cutoff Rule 
            cutoff = [nanmean(sVals)-nanstd(sVals) nanmean(sVals)+nanstd(sVals)];
            cutoff = [-inf inf];

            % Tukey plot rule +/- 1.6 IQR
%             iqr = [sVals(ceil(length(sVals).*0.75))-sVals(floor(length(sVals).*0.25))];
%             cutoff = [sVals(floor(length(sVals).*0.25))-iqr.*1.5 sVals(ceil(length(sVals).*0.75))+iqr.*1.5];

            plot(mean([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i))]).*ones(1,2),...
                [max(cutoff(1),sVals(1)) min(cutoff(2),sVals(end))],'color','k','linewidth',1.5,'linestyle','-');
            
%             plot(mean([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i))]),...
%                 sVals([sVals<cutoff(1)|sVals>cutoff(2)]),'linestyle','none',...
%                 'marker','d','markersize',2.5,'markerfacecolor','k','color','k','markeredgecolor','k')
            
            h(i,j) = patch([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i)) (j-(w./2))+(wpg.*(i)) (j-(w./2))+(wpg.*(i-1))],...
                [sVals(floor(length(sVals).*0.25)) sVals(floor(length(sVals).*0.25)) ...
                sVals(ceil(length(sVals).*0.75)) sVals(ceil(length(sVals).*0.75))],...
                groupColor(ceil(i),:),'edgecolor',groupColor(ceil(i),:),'linewidth',2);
            
            plot(([(j-(w./2))+(wpg.*(i-1)) (j-(w./2))+(wpg.*(i))]),...
                [nanmedian(sVals) nanmedian(sVals)],'color','k','linewidth',2);
        end
    end
end