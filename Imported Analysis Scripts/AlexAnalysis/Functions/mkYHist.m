function h = mkYHist(d)
    h = nan;
    
    w = 0.8;
    gap = 0.1;
    iw = [w - ((length(d(:,1))+1).*(gap.*w))]./(length(d(:,1)).*2);
    
    tmp = [(1-w)./2 (gap.*w)];
    for i = 1:length(d(:,1))
        tmp = [tmp iw.*2 (gap.*w)];
    end
    tmp = cumsum(tmp)-0.5;
    ix = nanmean([tmp(2:2:end-1)' tmp(3:2:end)'],2)';
    
    colors = [hsv(length(d(:,1))).*0.25] +0.25;
    
    if length(d(:,1))==1
        colors = [0.5 0.5 0.5];
    end

    ud = unique(cat(1,d{:}));
    ud = [0; ud];
    for j = 1:length(d(1,:))
        for i = 1:length(d(:,1))
            t = histc(d{i,j},ud);
            t = t./nanmax(t);
            h(i) = patch([j+ix(i)-t.*iw; j+ix(i)+flipud(t).*iw],[ud; flipud(ud)],colors(i,:),...
                'linestyle','none');
            hold on
        end
    end
    set(gca,'xlim',[0 length(d(1,:))+1],'xtick',1:length(d(1,:)));
end