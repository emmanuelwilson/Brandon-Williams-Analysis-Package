function plotEBCs(paths)
    for p = paths'
        s = load(p{1});
        
%         doInclude = s.processed.splithalf.ebc.mrl.p < 0.05 & ...
%             s.processed.splithalf.ebc.shc.p < 0.05 & ...
%             s.processed.exclude.SFPs;
        
        doInclude = s.processed.splithalf.ebc.mrl.p < 0.05 & ...
            s.processed.ebc.shPFD < 45 & s.processed.ebc.shD < 3 & ...
            s.processed.exclude.SFPs;
        
        m = s.processed.ebc.whole(:,:,doInclude);

        doK = [5 5];
        for part = 0:floor(length(m)/prod(doK))

            figure(1)
            set(gcf,'position',[50 50 800 800],'color','k')
            for k = 1:prod(doK)
                if part.*prod(doK)+k > length(m(1,1,:))
                    break
                end
                subplot(doK(1),doK(2),k)
                showEBC(m(:,:,part.*prod(doK)+k))
                hold on
%                 plot([[d(k)-0.5].*sind(pfd(k))],[[d(k)-0.5].*cosd(pfd(k))],'marker','o', ...
%                     'markeredgecolor','w');
                colormap('jet')
                axis equal
                axis off
                drawnow
            end
            
            slashInds = find(ismember(p{1},'/'));
            outP = ['Plots/EBCMaps/' p{1}(slashInds+1:end-4) '_Partition_' num2str(part+1)];
            saveFig(gcf,outP,'tiff')
            saveFig(gcf,outP,'pdf')
            close all
            drawnow
        end
    end
end