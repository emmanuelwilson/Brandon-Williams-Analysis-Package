function mkSessionComps(sim,envs,doComps,root)
    checkP(root);
    toPlot = repmat({[]},[1 length(doComps(:,1))]);
    lags = abs(bsxfun(@minus,[1:length(sim(:,1,1))]',[1:length(sim(:,1,1))]));
    compLabel = [];
    for i = 1:length(doComps(:,1))
        compLabel = [compLabel {[doComps{i,1} '-' doComps{i,2}]}];
        ga = bsxfun(@times,ismember(envs,doComps(i,1)),ismember(envs,doComps(i,2))');
        gb = bsxfun(@times,ismember(envs,doComps(i,2)),ismember(envs,doComps(i,1))');
        isGood = (ga|gb);
        if ~any(isGood)
            continue
        end
        inc = cat(1,sim(repmat(isGood,[1 length(sim(1,1,:))])));
        toPlot{i} = inc(~isnan(inc));
    end
    figure
    set(gcf,'position',[50 50 25.*length(doComps(:,1)) 350])
    mkGraph([toPlot(1:2:end); toPlot(2:2:end)],doComps(1:2:end,2)')
%     mkWhisker([toPlot(1:2:end); toPlot(2:2:end)],doComps(1:2:end,2)')
    set(gca,'ylim',[-1 1])
    hold on
    plot(get(gca,'xlim'),[0 0],'linewidth',1,'linestyle','--','color','k')
    xlabel('Condition Comparison')
    saveFig(gcf,root,[{'pdf'} {'tiff'}]);
end