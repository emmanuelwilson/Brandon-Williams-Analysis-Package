function addEBCmaps(paths,varargin)
    clc
    close all
    drawnow
    
    params = ebcMapParams();
    
    if isempty(varargin)
        varargin = [{'ebc'}];
    end
    varargin = cellfun(@lower,varargin,'uniformoutput',false);
    
    warning off all
    if isempty(gcp)
        parpool('local',7);
    end
    pctRunOnAll warning off all
    
    
    fprintf(['Computing EBC maps:\n']);
    tmp = [repmat({'\n\t'},[1 length(varargin)]); varargin];
    fprintf(cat(2,tmp{:},'\n\n'))
    for p = paths'
        s = load(p{1});

        
        fprintf(['\n\t' num2str(p{1})])    
        
%         if isfield(s.processed.ebc,'walls')
% %             continue
%         end
        
        fprintf('\t\tComputing EBC maps...  ')
        tic

        s.processed.p = interpNaNs(s.processed.p')';
        [m wm] = mkEBCMap(s.processed.p,s.processed.hd,s.processed.trace);

        s.processed.ebc.whole = m;
        s.processed.ebc.walls = wm;
        
        [mrl mrld pfd] = ebcScore(m);
        
        
        alpha = [2.*pi./(2.*length(m(:,1,1))):2.*pi./length(m(:,1,1)):2.*pi];
        
%         anulus = nan(length(m(:,1,1,1)),length(m(1,1,:,1)));
%         wanulus = nan(length(m(:,1,1,1)),length(m(1,1,:,1)),4);
%         for k = 1:length(m(1,1,:,1))
%             anulus(:,k) = m(:,mrld(k),k);
%             wanulus(:,k,:) = permute(wm(:,mrld(k),k,:),[1 2 4 3]);
%         end
%         pfd = rad2deg(circ_mean(repmat(alpha',[1 length(anulus(1,:))]),anulus));

        s.processed.ebc.mrl = mrl;
        s.processed.ebc.mrld = (mrld-1).*params.binDist + params.binDist./2;
        s.processed.ebc.pfd = rad2deg(alpha(pfd));
        
        wmc = nan(length(wm(1,1,1,:)),length(wm(1,1,1,:)),length(wm(1,1,:,1)));
        for i = 1:length(wm(1,1,1,:))
            for j = i+1:length(wm(1,1,1,:))
                t1 = wm(:,:,:,i);
                t2 = wm(:,:,:,j);
                t1 = reshape(t1,[numel(wm(:,:,1,1)) length(wm(1,1,:,1))]);
                t2 = reshape(t2,[numel(wm(:,:,1,1)) length(wm(1,1,:,1))]);
                inc = ~any(isnan(t1),2) & ~any(isnan(t2),2);
                xc = corr(t1(inc,:),t2(inc,:));
                wmc(i,j,:) = xc(logical(eye(size(xc))));
            end
        end
        
        s.processed.ebc.medianWallCorr =  nanmedian(reshape(wmc,[numel(wmc(:,:,1)) length(wmc(1,1,:))]),1)';
        
        mask = ([1:length(s.processed.hd)]<length(s.processed.hd)./2);
        m1 = mkEBCMap2(s.processed.p,s.processed.hd,s.processed.trace,{mask});
        m2 = mkEBCMap2(s.processed.p,s.processed.hd,s.processed.trace,{~mask});
        
        rm1 = reshape(m1,[numel(m1(:,:,1)) length(m1(1,1,:))]);
        rm2 = reshape(m2,[numel(m2(:,:,1)) length(m2(1,1,:))]);
        isGood = ~any(isnan(rm1),2) & ~any(isnan(rm2),2);
        xc = corr(rm1(isGood,:),rm2(isGood,:));
        xc = xc(logical(eye(size(xc))));
        s.processed.ebc.split = cat(4,m1,m2);
        s.processed.ebc.shc = xc;
        
        [a1 b1 c1] = ebcScore(m1);
        [a2 b2 c2] = ebcScore(m2);
        
        [a b] = getAngDiff(m1,m2);
        
%         tmp = mod([c1-c2].*(360./length(m1(:,1,1))),360);
        s.processed.ebc.shPFD = a; %nanmin(tmp,360-tmp);
        s.processed.ebc.shD = abs(b1-b2);
        
        durat = toc;
        fprintf([num2str(durat) ' sec']);
        
        save(p{1},'-struct','s','-v7.3');
    end
end

function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end

function sh = help_splithalf(p,t)

    v = [0 sqrt(sum(diff(p,[],2).^2))].*30;
    m1 = mkTraceMaps(p,t,1:length(p(1,:))<length(p(1,:))./2);
    m2 = mkTraceMaps(p,t,1:length(p(1,:))>=length(p(1,:))./2);

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