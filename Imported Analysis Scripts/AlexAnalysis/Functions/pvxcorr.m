function am = pvxcorr(m1,m2,lims,minSamples)
    if nargin==1 || isempty(m2)
        m2 = m1;
    end

    if any(size(m2)<size(m1))
        tmp = m2;
        m2 = m1;
        m1 = tmp;
    end

    if any(size(m2)~=size(m1))
        biggest = max(size(m2),size(m1));
        m1 = padarray(m1,floor([biggest-size(m1)]./2),nan,'both');
        m2 = padarray(m2,floor([biggest-size(m2)]./2),nan,'both');
        
        biggest = max(size(m2),size(m1));
        m1 = padarray(m1,([biggest-size(m1)]),nan,'pre');
        m2 = padarray(m2,([biggest-size(m2)]),nan,'pre');
    end
    
    if nargin<4
        minSamples = 0;
    end
    
    if nargin<3 || isempty(lims)
        lims = size(m2(:,:,1))-1;
    elseif all(lims<1)
        lims = round((size(m2(:,:,1))-1).*lims);
    end
    
    if length(lims)==1
        lims = [lims lims];
    end
    
    [a b c] = size(m1);
    am = nan(2.*a-1,2.*b-1);
    maxLims = (size(m2(:,:,1))-1);
    lims(lims>maxLims)=maxLims(lims>maxLims);
    for xlag = ceil(length(am(:,1,1))./2)-(lims(1)):ceil(length(am(:,1,1))./2)+(lims(1))
        xa = [max(length(m1(:,1,1))-xlag+1,1):...
            min(length(m1(:,1,1)),2.*length(m1(:,1,1))-xlag)];
        xb = [max(-length(m1(:,1,1))+xlag+1,1):...
            min(length(m1(:,1,1)),xlag)];
        for ylag = ceil(length(am(1,:,1))./2)-(lims(2)):ceil(length(am(1,:,1))./2)+(lims(2))
            ya = [max(length(m1(1,:,1))-ylag+1,1):...
                min(length(m1(1,:,1)),2.*length(m1(1,:,1))-ylag)];
            yb = [max(-length(m1(1,:,1))+ylag+1,1):...
                min(length(m1(1,:,1)),ylag)];
            ind1 = false(size(m1(:,:,1)));
            ind2 = false(size(m1(:,:,1)));
            ind1(xa,ya) = true;
            ind2(xb,yb) = true;
            
            vals = [m1(repmat(ind1,[1 1 length(m1(1,1,:))])),...
                m2(repmat(ind2,[1 1 length(m1(1,1,:))]))];
            vals(any(isnan(vals),2),:) = [];
            if isempty(vals)
                continue
            end
            if length(vals(:,1)) < minSamples
                continue
            end
            am(xlag,ylag) = corr(vals(:,1),vals(:,2));
        end
    end
    am = am(max(ceil(length(am(:,1,1))./2)-(lims(1)),1):min(ceil(length(am(:,1,1))./2)+(lims(1)),length(am(:,1,1))),...
        [max(ceil(length(am(1,:,1))./2)-(lims(2)),1):min(ceil(length(am(1,:,1))./2)+(lims(2)),length(am(1,:,1)))]);
end