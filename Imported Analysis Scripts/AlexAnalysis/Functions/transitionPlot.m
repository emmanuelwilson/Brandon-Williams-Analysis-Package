function h = transitionPlot(data,envs,doComps,doDots)
    if nargin <4 | isempty(doDots)
        doDots = true;
    end

    if ~iscell(data)
        data = num2cell(data);
    end
    data(logical(eye(size(data)))) = {nan};

    toPlot = repmat({[]},[1 length(doComps(:,1))]);
    compLabel = [];
    for i = 1:length(doComps(:,1))
        compLabel = [compLabel {[doComps{i,1} '-' doComps{i,2}]}];
        ga = bsxfun(@times,ismember(envs,doComps(i,1)),ismember(envs,doComps(i,2))');
        gb = bsxfun(@times,ismember(envs,doComps(i,2)),ismember(envs,doComps(i,1))');
        isGood = (ga|gb);
        if ~any(isGood)
            continue
        end
        inc = cat(1,data{isGood});
        toPlot{i} = inc(~isnan(inc));
    end
    pC = [];
    for i = 1:length(doComps(:,1))
        pC = [pC; {[doComps{i,1} '-' doComps{i,2}]}];
    end
%     figure
    set(gcf,'position',[50 50 60.*length(doComps(:,1)) 350])
    h = mkGraph([toPlot],pC',doDots);
%     set(gca,'ylim',[-0.5 1])
    hold on
    plot(get(gca,'xlim'),[0 0],'linewidth',1,'linestyle','-','color','k')
    xlabel('Condition Comparison')
end