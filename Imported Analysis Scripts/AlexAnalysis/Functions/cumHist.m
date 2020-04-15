function [h at] = cumHist(d,x)
    h = [];

    if isnumeric(d)
        tmp = [];
        for i = 1:length(d(1,:))
            tmp = [tmp {d(:,i)}];
        end
        d = tmp;
    end

    set(gca,'xlim',[x(1) x(end)],'ylim',[0 1])
%     plot([0 1],[0 1],'linestyle','--','linewidth',2,'color',[0.6 0.6 0.6])
    hold on
%     x = [0:0.025:1];
%     colors = transcm;
%     colors = colors(([1 40 160 200]),:);
%     colors = [0.3 0.3 0.9; 0.7 .7 .9; 0.9 0.3 0.3; 0.9 0.7 0.7];
    colors = [0.0 0.0 0.0; 0.3 0.3 0.3; 0.5 0.5 0.5; 0.7 0.7 0.7; 0.8 0.8 0.8; 0.3 0.3 0.3];
%     colors = hsv(length(d));

%     colors = repmat([0.5 0.5 0.5],[length(d) 1]);
%     colors = repmat([1 0.2 0.2],[length(d) 1]);
    
%     colors = ([0.9 0.1 0.1; 0.65 0.1 0.1; 0.5 0.1 0.1; 0.35 0.1 0.1; 0.2 0.1 0.1; 0.05 0.05 0.05]);

%     colors = [0.25 0.25 0.25; 0.0 0.9 1; 0.25 0.25 0.25; 0.0 0.9 1; 0.25 0.25 0.25; 0.0 0.9 1 ; ...
%         0.25 0.25 0.25; 0.0 0.9 1; 0.25 0.25 0.25; 0.0 0.9 1; 0.25 0.25 0.25; 0.0 0.9 1];

    at = [];
    ls = [{'-'} {'--'} {'--'} {'--'} {'-'} {'--'} {'-'} {'-'} {'-'} {'-'} {'-'} {'-'}];
    for i = 1:length(d)
        if isempty(d{i})
            continue
        end
        t = histc(d{i},x);
        t(end) = [];
        t = cumsum(t);
        t = t./(t(end));
        at = [at t];
        h{i} = plot([x(1) x(1:end-1)+((x(2)-x(1))/2) x(end)],[0 t(:)' 1],'linewidth',2,...
            'color',colors(i,:).*0.85,'linestyle',ls{i});
        hold on
    end
end