function similarityAcrossSessions(folder)
    mice = dir(folder);
    mice = {mice(3:end).name};
    
    similarity = repmat({[]},[1 length(mice)]);
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
        
        s1 = load([mp '/' sessions{1}]);
        similarity{ismember(mice,m)} = nan(length(sessions),length(sessions),nansum(s1.processed.isAligned));
        for dci = 1:length(sessions)
            s1 = load([mp '/' sessions{dci}]);
            for dcj = dci+1:length(sessions)
                s2 = load([mp '/' sessions{dcj}]);

                v = [0 sqrt(sum(diff(s1.processed.p,[],2).^2))].*30;
                m1 = mkTraceMaps(s1.processed.p,s1.processed.trace(s1.processed.isAligned,:),[]);
                v = [0 sqrt(sum(diff(s2.processed.p,[],2).^2))].*30;
                m2 = mkTraceMaps(s2.processed.p,s2.processed.trace(s2.processed.isAligned,:),[]);
                big = nanmax([size(m1); size(m2)]);

                v = [0 sqrt(sum(diff(s1.processed.p,[],2).^2))].*30;
                m1 = mkTraceMaps(s1.processed.p,s1.processed.trace(s1.processed.isAligned,:),[],big(1:2));
                v = [0 sqrt(sum(diff(s2.processed.p,[],2).^2))].*30;
                m2 = mkTraceMaps(s2.processed.p,s2.processed.trace(s2.processed.isAligned,:),[],big(1:2));

                actual = m2sh(m1,m2);   
                similarity{ismember(mice,m)}(dci,dcj,:) = actual;
%                 nsims = 1000;
%                 null = nan(length(s1.processed.trace(s1.processed.isAligned,1)),nsims);
%                 parfor sim = 1:nsims
%                     null(:,sim) = m2sh(m1,m2(:,:,randperm(length(m2(1,1,:)))));   
%                 end
% 
%                 pval = 1-nanmean(bsxfun(@gt,actual',null(:)))';
            
            end
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