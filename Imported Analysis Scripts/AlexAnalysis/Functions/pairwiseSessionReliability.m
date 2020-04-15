
function pairwiseSessionReliability()
    
    clc
    fprintf('\n')
    load('combinedData_OpenField');
     
%     dat(3) = [];
    %%% Reliability constaint

    for mi = 1:length(dat)
        dat(mi).rely = repmat({[]},[6 6]);
        dat(mi).doskip = false(6,6);
        fprintf(['\n\tMouse:  ' num2str(mi) '\n'])
        
        for di = 1:1:length(dat(mi).maps.overall(1,:))
            for dj = di+1:1:length(dat(mi).maps.overall(1,:))
            
                %%% Test all rotations
                vals = nan(1,4);
                for rot = [0 2]
                    doM1 = dat(mi).maps.overall{:,di};
                    doM2 = dat(mi).maps.overall{:,dj};
                    
                    doM1 = doM1(:,:,dat(mi).registration{di,dj}(:,1));
                    doM2 = doM2(:,:,dat(mi).registration{di,dj}(:,2));
                    
                    doM2 = rot90(doM2,rot);

                    isGood = ~isnan(doM1)&~isnan(doM2);
                    if ~any(isGood(:))
                        continue
                    end
                    vals(rot+1) = corr(doM1(isGood),doM2(isGood));
                end

                [blah bestRot] = nanmax(vals);

                doM1 = dat(mi).maps.overall{:,di};
                doM2 = dat(mi).maps.overall{:,dj};

                doM1 = doM1(:,:,dat(mi).registration{di,dj}(:,1));
                doM2 = doM2(:,:,dat(mi).registration{di,dj}(:,2));
                    
                doM2 = rot90(doM2,bestRot-1);
                reliability = nan(length(doM1(1,1,:)),1);
                for k = 1:length(doM1(1,1,:))
                    m1 = doM1(:,:,k);
                    m2 = doM2(:,:,k);
                    isGood = ~isnan(m1)&~isnan(m2);
                    reliability(k) = corr(m1(isGood),m2(isGood));
                end

                nsims = 250;
                null = nan(length(reliability),nsims);
                for si = 1:nsims
                    doM1 = doM1(:,:,randperm(length(doM1(1,1,:))));
                    for k = 1:length(doM1(1,1,:))
                        m1 = doM1(:,:,k);
                        m2 = doM2(:,:,k);
                        isGood = ~isnan(m1)&~isnan(m2);
                        null(k,si) = corr(m1(isGood),m2(isGood));
                    end
                end

                [a b] = sort(null(:));
                thresh = a(round(length(a).*0.95));
                dat(mi).rely{di,dj} = reliability >= thresh;
                dat(mi).rely{di,dj} = reliability >= thresh;
                dat(mi).rotation{di,dj} = bestRot-1;

                fprintf(['\t\t( ' num2str(di) ', ' num2str(dj) '; rot:  ' num2str((bestRot-1).*90) ') Count:  ' num2str(nansum(reliability >= thresh)) ...
                    '\tProportion:  ' num2str(nanmean(reliability >= thresh)) '\n'])

                if nansum(reliability >= thresh) < 40 %%% Toss sessions with fewer than 30 reliable place cells
                    dat(mi).doskip(di,dj) = true;
                    dat(mi).doskip(di,dj) = true;
                else
                    dat(mi).doskip(di,dj) = false;
                    dat(mi).doskip(di,dj) = false;
                end
            end
        end
    end
end


function vals = help_getMaskedVals(v,mask)
    vals = nan(length(v(1,1,:)),1);
    for k = 1:length(v(1,1,:))
        tmp = v(:,:,k);
        vals(k) = nanmean(tmp(mask));
    end
end















