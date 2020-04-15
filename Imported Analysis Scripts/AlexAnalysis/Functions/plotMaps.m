function plotMaps(paths)
    for p = paths'
        s = load(p{1});
        
        v = [0 sqrt(sum(diff(s.processed.p,[],2).^2))].*30;
        m = mkTraceMaps(s.processed.p,s.processed.trace,[]);
        
        m(:,:,s.processed.splithalf.wholemap_unmatched.p > 0.05) = [];
        
        doK = [8 8];
        for part = 0:floor(length(m)/prod(doK))

            figure(1)
            set(gcf,'position',[50 50 900 900],'color','k')
            for k = 1:prod(doK)
                if part.*prod(doK)+k > length(m(1,1,:))
                    break
                end
                subplot(doK(1),doK(2),k)
                imagesc(m(:,:,part.*prod(doK)+k))
                colormap('jet')
                alpha(double(~isnan(m(:,:,part.*prod(doK)+k))))
%                 title(num2str(nanmax(nanmax(m(:,:,part.*prod(doK)+k)))))
%                 drawnow
%                 set(gca,'ydir','normal')
                axis equal
                axis off
            end
            
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/TraceFieldMaps/' p{1}(slashInds+1:end-4) '_Partition_' num2str(part+1)];
            saveFig(gcf,outP,'tiff')
            saveFig(gcf,outP,'pdf')
            close all
            drawnow
        end
    end
end