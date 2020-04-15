function plotMultiNDayMaps(paths)
    
    warning off all
    
    clc
    fprintf('\n')
    %%% Reliability constaint
    
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
    
    doPlot = true;
    
    envLabel = [{'square'} {'glenn'}];
    for mi = 1:length(upiece)
        fprintf(['\n\tMouse:  ' num2str(upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end)) '\n'])
        isM = find(ismember(piece,upiece(mi)));
        sessions = paths(isM);
        s = load(paths{isM(1)},'alignment');
        
        %%% Check to confirm aligned paths are plotted paths
        
        if ~all(cellfun(@strcmp,s.alignment.sessions,sessions))
            fprintf('\n\t****** ERROR - ALIGNED SESSIONS ARE NOT PLOTTED SESSIONS ******');
        end
        
        alignMap = s.alignment.alignmentMap;
        combs = s.alignment.combs;
        am = repmat({[]},[1 length(sessions)]);
        isPC = repmat({[]},[1 length(sessions)]);
        
        clear s
        for si = 1:length(sessions)
            fprintf(['\t\tSession:  ' sessions{si}(find(ismember(sessions{si},'/'),1,'last')+1:end) '\n'])
            s(si) = load(sessions{si},'processed','properties');
        end 
        
        for ci = 1:length(combs(:,1))
            fprintf(['\t\t\tCombination:  ' num2str(combs(ci,:)) '\n'])
            cP = [];
            cT = [];
            env = [];
            days = [];
            for cj = 1:length(combs(1,:))
                cP = [cP s(combs(ci,cj)).processed.p];
                cT = [cT s(combs(ci,cj)).processed.trace(alignMap{ci}(:,cj),:)];
                cEnv = lower(sessions{combs(ci,cj)}(find(ismember(sessions{combs(ci,cj)},'_'),1,'last')+1:end-4));
                env = [env repmat(find(ismember(envLabel,cEnv))-1,[1 length(s(combs(ci,cj)).processed.p(1,:))])];
                days = [days repmat(cj,[1 length(s(combs(ci,cj)).processed.p(1,:))])];
            end
            allMasks = repmat({[]},[1 nanmax(days)]);
            for i = 1:nanmax(days)
                allMasks{i} = days==i;
            end
            
            clear map
            [map samp allComp] = getMatchedMapsNMasks(cP,cT,allMasks);
            
            map = correctRotation(map);
            
            if doPlot
                outP = ['Plots/NWiseAlignedCellMaps/' upiece{mi}(find(ismember(upiece{mi},'/'),1,'last')+1:end) '/' ...
                    'Combination_' num2str(ci)];
                helpPlotMaps(map,outP);
            end
        end
    end
end

function helpPlotMaps(m,outP)
    doK = [8 4];

    toPlot = [];
    for i = 1:length(m(1,1,1,:))
        if i > 1
            toPlot = cat(2,toPlot,nan(size(m(:,1,:,1))));
        end
        norm = repmat(nanmax(nanmax(m(:,:,:,i),[],1),[],2),[size(m(:,:,1,i))]);
        toPlot = cat(2,toPlot,m(:,:,:,i)./norm);
    end
    
    for part = 0:floor(length(toPlot(1,1,:,1))/prod(doK))

        figure(1)
        set(gcf,'position',[50 50 900 1350])
        for k = 1:prod(doK)
            if part.*prod(doK)+k > length(toPlot(1,1,:))
                break
            end

            subplot(doK(1),doK(2),k)
            imagesc(toPlot(:,:,part.*prod(doK)+k))
            colormap jet
            alpha(double(~isnan(toPlot(:,:,part.*prod(doK)+k))))
            axis equal
            axis off    
        end
        outP2 = [outP '_Part_' num2str(part)];
        saveFig(gcf,outP2,[{'tiff'} {'pdf'}])
        close all
        drawnow
    end
end