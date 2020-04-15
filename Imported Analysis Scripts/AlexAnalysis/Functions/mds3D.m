function mds3D(mat,envs,envLabel,root)

    %%% Multidimensional 3D video
%     tmp = cellfun(@nanmedian,kCrossSim);
    if iscell(mat)
        tmp = cellfun(@nanmedian,mat);
    else
        tmp = nanmedian(mat,3);
    end
    tmp = nanmax(tmp,tmp');
    tmp(logical(eye(size(tmp)))) = 1;
    tmp = 1-tmp;

    [mdssim eigenVals] = cmdscale(tmp);
    
%     eigenVals./nanmax(eigenVals);
%     fprintf(['\n\tMax relative error (1D)= ' num2str(nanmean(abs(tmp - squareform(pdist(mdssim(:,1))))) / max(tmp))]);
%     fprintf(['\n\tMax relative error (2D)= ' num2str(nanmean(abs(tmp - squareform(pdist(mdssim(:,1:2))))) / max(tmp))]);
%     fprintf(['\n\tMax relative error (3D)= ' num2str(nanmean(abs(tmp - squareform(pdist(mdssim(:,1:3))))) / max(tmp))]);
%     fprintf(['\n\tMax relative error (4D)= ' num2str(nanmean(abs(tmp - squareform(pdist(mdssim(:,1:4))))) / max(tmp))]);
    mdssim = mdssim(:,1:3);
    
    figure
    set(gcf,'position',[50 50 800 800])
    lim = nanmax(abs(mdssim(:)))+0.1;
    set(gca,'xlim',[-lim lim],'ylim',[-lim lim],'zlim',[-lim lim])
    hold on
    keyShape = ['oooooo'];
    keyColor = [0.9 0.1 0.1; 0.9 0.4 0.4; 0.9 0.7 0.7; 0.7 0.7 0.9; 0.4 0.4 0.9; 0.1 0.1 0.9].*0.5+0.5;
    initSize = 6;
    keySize = initSize+[floor([1:length(mat)]./6)].*2.5;
    for i = 1:length(envLabel)
        envI = find(ismember(envs,envLabel(i)));
        plot3(mdssim(envI,1),mdssim(envI,2),mdssim(envI,3),'marker','none', ...
            'color',keyColor(i,:),'linewidth',1);
    end

    allH = [];
    for i = 1:length(mdssim(:,1))
        envI = find(ismember(envLabel,envs(i)));
        h = plot3(mdssim(i,1),mdssim(i,2),mdssim(i,3),'marker',keyShape(envI), ...
            'color',keyColor(envI,:),'markerfacecolor',keyColor(envI,:),'markersize',keySize(i),...
            'markeredgecolor','w');
        if i == nanmax(find(ismember(envs,envLabel(envI))))
            allH(envI) = h;
            set(h,'markeredgecolor','k','linewidth',1.5);
%                 text(mdssim(i,1),mdssim(i,2),upper(envLabel(envI)),'fontname','arial',...
%                     'fontsize',6,'fontweight','bold','horizontalalignment','center',...
%                     'verticalalignment','middle')
        end
    end


    legend(allH,[{'Square'} {''} {''} {''} {''} {'Glenn'}],...
        'location','none','orientation','horizontal','box','off','position',[0.1 -0.025 0.8 0.1])
    zlabel('MDS Dim 3');
    ylabel('MDS Dim 2');
    xlabel('MDS Dim 1');
%         axis equal
    grid on
    clear frames
    for i = 1:1:360
        view(i+0.5,20)
        drawnow
        frames(i) = getframe(gcf);

        im = frame2im(frames(i)); 
        [imind,cm] = rgb2ind(im,256); 
        if i == 1
            checkP(root);
            imwrite(imind,cm,root,'gif', 'Loopcount',inf,'DelayTime',0.025); 
        else
            imwrite(imind,cm,root,'gif','DelayTime',0.025,'WriteMode','append'); 
        end
    end
end