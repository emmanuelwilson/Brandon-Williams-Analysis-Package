function [] = SectionDisplay(ms,frameMap,times)
mkdir SectionCalcium
t1 =frameMap(times(:,1));
t2 =frameMap(times(:,2));
times = cat(2,t1,t2);
traces = ms.FiltTraces';
for i = 1:length(times(:,1))
    for j = 1 : round(length(ms.FiltTraces(1,:)))/50
        figure(1)
        if (j == round(length(ms.FiltTraces(1,:)))/50)
            X= ms.FiltTraces(times(i,1):times(i,2),(j*50-50+1):end)'; % chose Ca traces
        else
            X= ms.FiltTraces(times(i,1):times(i,2),(j*50-50+1):j*50)'; % chose Ca traces
        end
        z = 0.11; % spacing parameter - experiment with different values
        vpos = z + z.*[1:size(X,1)]'; % calculates the spacing
        X1 = X + repmat(vpos,1,size(X,2)); % adds the spacing to the matrix
        l = size(X1,2); %Define number of columns
        m = size(X1,1); %Define number of rows
        fig = figure('color','w');
        hold on;
        set(fig,'DefaultLineLineWidth',0.8); %Ca trace linewidth
        jplot(X1);
        set(gca,'FontName','Arial','Fontsize', 14);
        xlabel('Time (s)')
        ylabel('Cell number');
        set(gca, 'box', 'off');
        set(gca,'TickDir','out');
        ylim([0 (size(X,1)*z+2*z)]);
        %xlim([0 8000]);
        set(gca,'XTick',[0 l*.25 l*.5 l*.75 l],'XTickLabel',{'0', ((l/30)*.25), ((l/30)*.5), ((l/30)*.75), l/30});
        set(gca,'YTick',[2*z m*z+z],'YTickLabel',{'1', m});
        set(gca, 'LineWidth', 1);
        
        saveas(gcf,['SectionCalcium','/',num2str(times(i,1)),'-',num2str(times(i,2)),'Cells',num2str(j*50),'.jpg'])
    end
end
end


