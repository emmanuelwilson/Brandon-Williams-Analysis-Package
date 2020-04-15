function saveFig(h,p,form)
    
    if ~iscell(form)
        form = {form};
    end

    checkP(p);
    
    ax = findall(h,'type','axes');
    set(ax,'fontname','arial','fontweight','bold','fontsize',11);
    drawnow;
    
    checkP(p);
    STYLE = hgexport('readstyle','default');
    
    for i = 1:length(form)
        hgexport(h, [p '.' form{i}],STYLE, 'Format',form{i});
    end
end