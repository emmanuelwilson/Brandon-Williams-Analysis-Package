function mds2D_binned(mat,envs,envLabel,numBins,root)
    %%% Multidimensional Scaling 2D
    if iscell(mat)
        tmp = cellfun(@nanmedian,mat);
    else
        tmp = nanmedian(mat,3);
    end
    tmp = nanmax(tmp,tmp');
    tmp(logical(eye(size(tmp)))) = 1;
    tmp = 2-(tmp).*2;

    [mdssim eigenVals] = cmdscale(tmp);
    mdssim = mdssim(:,1:2);

    features = [];
    for i = 1:length(envLabel)
        envI = find(ismember(envs,envLabel(i)));
        for k = 1:length(envI)
            features((envI(k)-1).*numBins+1:(envI(k)).*numBins,:) = ...
                [repmat([k i],[numBins 1]) [1:numBins]'];
        end
    end

    figure
    set(gcf,'position',[50 50 400 400])
    lim = nanmax(abs(mdssim(:)))+0.1;
    set(gca,'xlim',[-lim lim],'ylim',[-lim lim])
    hold on
    keyShape = ['oooooo'];
    keyColor = [0.9 0.1 0.1; 0.9 0.4 0.4; 0.9 0.7 0.7; 0.7 0.7 0.9; 0.4 0.4 0.9; 0.1 0.1 0.9].*0.5+0.5;
    initSize = 6;
    keySize = initSize+[floor([1:length(tmp)]./(6.*numBins))].*1;
    for i = 1:length(envLabel)
        envI = find(features(:,2)==i);
        plot(mdssim(envI,1),mdssim(envI,2),'marker','none', ...
            'color',keyColor(i,:),'linewidth',1);
    end

    for i = 1:length(mdssim(:,1))
        envI = features(i,2);
        h = plot(mdssim(i,1),mdssim(i,2),'marker',keyShape(envI), ...
            'color',keyColor(envI,:),'markerfacecolor',keyColor(envI,:),'markersize',keySize(i),...
            'markeredgecolor','w');
        if i == nanmax(find(features(:,2)==envI))
            set(h,'markeredgecolor','k','linewidth',1.5);
%             text(mdssim(i,1),mdssim(i,2),upper(envLabel(envI)),'fontname','arial',...
%                 'fontsize',6,'fontweight','bold','horizontalalignment','center',...
%                 'verticalalignment','middle')
        end
    end
    ylabel('MDS Dim 2');
    xlabel('MDS Dim 1');
    axis square
    
    checkP(root);
    saveFig(gcf,[root],[{'pdf'} {'tiff'}]);
end