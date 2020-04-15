function selectionAcrossSessions(folder)
    mice = dir(folder);
    mice = {mice(3:end).name};
    doComparisons = [2 3];
    for m = mice
        mp = [folder '/' m{1}];
        sessions = dir(mp);
        sessions = {sessions(3:end).name};
        
        % order sessions by date
        clear dates
        for j = 1:length(sessions)
            dates(j) = datetime(str2num(sessions{j}(1:2)),str2num(sessions{j}(4:5)),str2num(sessions{j}(7:8)));
        end
        [a b] = sort(dates);
        sessions = sessions(b);
        
        s1 = load([mp '/' sessions{doComparisons(1)}]);
        s2 = load([mp '/' sessions{doComparisons(2)}]);
        
%         v = [0 sqrt(sum(diff(s1.processed.p,[],2).^2))].*30;
        m1 = mkTraceMaps(s1.processed.p,s1.processed.trace,[]);
%         v = [0 sqrt(sum(diff(s2.processed.p,[],2).^2))].*30;
        m2 = mkTraceMaps(s2.processed.p,s2.processed.trace,[]);
        big = nanmax([size(m1); size(m2)]);

%         v = [0 sqrt(sum(diff(s1.processed.p,[],2).^2))].*30;
        m1 = mkTraceMaps(s1.processed.p,s1.processed.trace,[],big(1:2));
%         v = [0 sqrt(sum(diff(s2.processed.p,[],2).^2))].*30;
        m2 = mkTraceMaps(s2.processed.p,s2.processed.trace,[],big(1:2));
        
        actual = m2sh(m1,m2);   
        
        nsims = 1000;
        null = nan(length(s1.processed.trace(:,1)),nsims);
        parfor sim = 1:nsims
            null(:,sim) = m2sh(m1,m2(:,:,randperm(length(m2(1,1,:)))));   
        end
        
        pval = 1-nanmean(bsxfun(@gt,actual',null(:)))';
        
        
%         actual = help_splithalf(s1.processed.p,s1.processed.trace,...
%             s2.processed.p,s2.processed.trace,big(1:2));
%       
%         nsims = 100;
%         null = nan(length(s1.processed.trace(:,1)),nsims);
%         parfor sim = 1:nsims
%             sp1 = circshift(s1.processed.p,[0 900+randi(length(s1.processed.p(1,:))-1800)]);
%             sp2 = circshift(s2.processed.p,[0 900+randi(length(s2.processed.p(1,:))-1800)]);
%             null(:,sim) = help_splithalf(sp1,s1.processed.trace,...
%             sp2,s2.processed.trace,big(1:2));
%         end
%         
%         pval= 1-nanmean(repmat(actual',[nsims 1]) > null')';
        
        figure(1)
        set(gcf,'position',[50 50 250 250])
        cumHist(pval,[0:0.01:1])
        hold on
        plot([0 1],[0 1],'color',[0.5 0.5 0.5],'linestyle',':','linewidth',1)
        ylabel('Cumulative Proportion')
        xlabel('P-Value')
        drawnow
        slashInds = find(ismember(mp,'/'));
        outP = ['Plots/StabilityAcrossSessions/' mp(slashInds(end)+1:end)];
        saveFig(gcf,outP,'tiff')
        close all
        drawnow

        
        %%% Save successfully registered cells
        for j = 1:length(sessions)
            ref = load([mp '/' sessions{j}]);
            
            
            ref.processed.splithalf.val = actual;
            ref.processed.splithalf.null = null;
            ref.processed.splithalf.p = pval;
            
            save([mp '/' sessions{j}],'-struct','ref','-v7.3');
        end
    end
end

function sh = m2sh(m1,m2)
    rm1 = nan(length(m1(1,:,1)).*length(m1(:,1,1)),length(m1(1,1,:)));
    rm2 = nan(length(m1(1,:,1)).*length(m1(:,1,1)),length(m1(1,1,:)));
    for k = 1:length(m1(1,1,:))
        tmp = m1(:,:,k);
        rm1(:,k) = tmp(:);
        tmp = m2(:,:,k);
        rm2(:,k) = tmp(:);
    end

    isBad = isnan(rm1) | isnan(rm2);

    m1t = reshape(rm1(~isBad),[nansum(~isBad(:,1)) length(m1(1,1,:))]);
    m2t = reshape(rm2(~isBad),[nansum(~isBad(:,1)) length(m1(1,1,:))]);

    xc = (corr(m1t,m2t));
    sh = xc(logical(eye(size(xc))));
end

function sh = help_splithalf(p1,t1,p2,t2,big)

%     v = [0 sqrt(sum(diff(p1,[],2).^2))].*30;
    m1 = mkTraceMaps(p1,t1,[],big);
%     v = [0 sqrt(sum(diff(p2,[],2).^2))].*30;
    m2 = mkTraceMaps(p2,t2,[],big);

    rm1 = nan(length(m1(1,:,1)).*length(m1(:,1,1)),length(m1(1,1,:)));
    rm2 = nan(length(m1(1,:,1)).*length(m1(:,1,1)),length(m1(1,1,:)));
    for k = 1:length(m1(1,1,:))
        tmp = m1(:,:,k);
        rm1(:,k) = tmp(:);
        tmp = m2(:,:,k);
        rm2(:,k) = tmp(:);
    end

    isBad = isnan(rm1) | isnan(rm2);

    m1t = reshape(rm1(~isBad),[nansum(~isBad(:,1)) length(m1(1,1,:))]);
    m2t = reshape(rm2(~isBad),[nansum(~isBad(:,1)) length(m1(1,1,:))]);

    xc = (corr(m1t,m2t));
    sh = xc(logical(eye(size(xc))));
end