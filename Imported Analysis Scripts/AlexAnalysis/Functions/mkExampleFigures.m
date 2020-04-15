function mkExampleFigures(paths)
    clc
    close all
    drawnow
    fprintf('\n')
    %%% Reliability constaint
    
    warning off all
    if isempty(gcp)
        parpool('local',7);
    end
    pctRunOnAll warning off all
    
    %% Split by animal
    piece = [];
    ag = [];
    spiece = [];
    for i = 1:length(paths)
        ind = find(ismember(paths{i},'/'),1,'last')-1;
        piece = [piece; {paths{i}(1:ind)}];
        spiece = [spiece; {paths{i}(ind+2:end-4)}];
    end
    upiece = unique(piece);
    
    for mi = 1:length(upiece)
        fprintf(['\n\n\tMouse:  ' num2str(upiece{mi}) '\n']) 
        isM = find(ismember(piece,upiece(mi)));
        for si = 1:length(isM);
            fprintf(['\n\t\tSession:  ' paths{isM(si)}])
            s = load(paths{isM(si)},'processed','properties','calcium');
            
            figure(1)
            set(gcf,'position',[50 50 800 200])
            subplot(1,4,1)
            imagesc(s.calcium.meanFrame(:,:,1))
            set(gca,'xlim',[0.5 251.5],'ylim',[0.5 160.5])
            colormap gray
            caxis([0 255])
            axis equal
            axis off
            if isfield(s.processed,'exclude')
                text(0,-20,sprintf(['Median SNR: %0.2f'], ...
                    nanmedian(s.processed.snr.whole(s.processed.exclude.SFPs))),...
                    'fontname','arial','fontweight','bold','fontsize',9)
            else
                text(0,-20,sprintf(['Median SNR: %0.2f'],nanmedian(s.processed.snr.whole)),...
                    'fontname','arial','fontweight','bold','fontsize',9)
            end
            text(0,180,['Mouse: ' upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end) '  ' ...
                'Session: ' s.properties.session],...
                'fontname','arial','fontweight','bold','fontsize',9)
            subplot(1,4,2)
            SFPs = s.calcium.SFPs;
            SFPs = SFPs./repmat(nanmax(nanmax(SFPs,[],1),[],2),[size(SFPs(:,:,1))]);
            if isfield(s.processed,'exclude')
                doSFPs = s.processed.exclude.SFPs;
            else
                doSFPs = true(length(SFPs(1,1,:)),1);
            end
            imagesc((nanmax(SFPs(:,:,doSFPs),[],3)).^(1.5))
            
            set(gca,'xlim',[0.5 251.5],'ylim',[0.5 160.5])
            axis equal
            axis off
            subplot(1,4,3:4)
            nK = 10;
            tmp = s.calcium.FiltTraces;
            if isfield(s.processed,'exclude')
                tmp = tmp(:,s.processed.exclude.SFPs);
            end
            doK = randperm(length(tmp(1,:)));
            plot(bsxfun(@plus,tmp(30.*60.*3:30.*60.*9,doK(1:nK)),[1:nK].*0.5), ...
                'linewidth',1,'color','k')
            hold on
            plot([-180 -180],[0.5 1.5],'color','k','linewidth',2)
            set(gca,'ylim',[-3 9],'xlim',[-210 30.*60.*6 + 60])
            axis off
            
            slashInds = find(ismember(paths{isM(si)},'/'));
            outP = ['Plots/TraceAndFrameExamples/' paths{isM(si)}(slashInds+1:end-4) ];
            saveFig(gcf,outP,[{'pdf'} {'tiff'}])
            close all
        end
    end
end