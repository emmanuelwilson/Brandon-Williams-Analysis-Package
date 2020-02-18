%% Compares EBCs Ratemaps B context EBCs across A context EBCs
function [out] = EBCstabilityratemapcomp_AABA_autocorrelation(Singlemap,fullpass1,fullpass2,fullpass3, fullpass4, ebc1, ebc2, ebc3, ebc4)
seshnum = 4;

[~, ind1, ~] = intersect(Singlemap(:,1),fullpass1);
[~, ind2, ~] = intersect(Singlemap(:,2),fullpass2);
[~, ind3, ~] = intersect(Singlemap(:,3),fullpass3);
[~, ind4, ~] = intersect(Singlemap(:,4),fullpass4);

bmap = Singlemap;
bmap(bmap>0) = 1;

CorrMap = cell(size(Singlemap,1),size(Singlemap,2));
maxCorrMap = nan(size(Singlemap,1),size(Singlemap,2));
medianCorrMap = nan(size(Singlemap,1),size(Singlemap,2));
meanCorrMap = nan(size(Singlemap,1),size(Singlemap,2));
minCorrMap = nan(size(Singlemap,1),size(Singlemap,2));
maxrotCorrMap = nan(size(Singlemap,1),size(Singlemap,2));
minrotCorrMap = nan(size(Singlemap,1),size(Singlemap,2));

for mapind = 1 : length(Singlemap(:,1))
    if sum(bmap(mapind,:)) > 1
        for sessind1 = 1 : seshnum
            if Singlemap(mapind,sessind1) > 0
                if sessind1 == 1
                    if Singlemap(mapind,sessind1) > 0
                        CorrMap{mapind,sessind1} = rmRotcorr_singlecell(ebc1.rm(:,:,Singlemap(mapind,sessind1)),ebc1.rm(:,:,Singlemap(mapind,sessind1)));
                        [maxCorrMap(mapind,sessind1), maxrotCorrMap(mapind,sessind1)] = nanmax(CorrMap{mapind,sessind1});
                        medianCorrMap(mapind,sessind1) = nanmedian(CorrMap{mapind,sessind1});
                        meanCorrMap(mapind,sessind1) = nanmean(CorrMap{mapind,sessind1});
                        [minCorrMap(mapind,sessind1), minrotCorrMap(mapind,sessind1)] = nanmin(CorrMap{mapind,sessind1});
                    end
                elseif sessind1 == 2
                    if Singlemap(mapind,sessind1) > 0
                        CorrMap{mapind,sessind1} = rmRotcorr_singlecell(ebc2.rm(:,:,Singlemap(mapind,sessind1)),ebc2.rm(:,:,Singlemap(mapind,sessind1)));
                        [maxCorrMap(mapind,sessind1), maxrotCorrMap(mapind,sessind1)] = nanmax(CorrMap{mapind,sessind1});
                        medianCorrMap(mapind,sessind1) = nanmedian(CorrMap{mapind,sessind1});
                        meanCorrMap(mapind,sessind1) = nanmean(CorrMap{mapind,sessind1});
                        [minCorrMap(mapind,sessind1), minrotCorrMap(mapind,sessind1)] = nanmin(CorrMap{mapind,sessind1});
                    end                    
                elseif sessind1 == 3
                    if Singlemap(mapind,sessind1) > 0
                        CorrMap{mapind,sessind1} = rmRotcorr_singlecell(ebc3.rm(1:25,:,Singlemap(mapind,sessind1)),ebc3.rm(1:25,:,Singlemap(mapind,sessind1)));
                        [maxCorrMap(mapind,sessind1), maxrotCorrMap(mapind,sessind1)] = nanmax(CorrMap{mapind,sessind1});
                        medianCorrMap(mapind,sessind1) = nanmedian(CorrMap{mapind,sessind1});
                        meanCorrMap(mapind,sessind1) = nanmean(CorrMap{mapind,sessind1});
                        [minCorrMap(mapind,sessind1), minrotCorrMap(mapind,sessind1)] = nanmin(CorrMap{mapind,sessind1});
                    end
                elseif sessind1 == 4
                    if Singlemap(mapind,sessind1) > 0
                        CorrMap{mapind,sessind1} = rmRotcorr_singlecell(ebc4.rm(:,:,Singlemap(mapind,sessind1)),ebc4.rm(:,:,Singlemap(mapind,sessind1)));
                        [maxCorrMap(mapind,sessind1), maxrotCorrMap(mapind,sessind1)] = nanmax(CorrMap{mapind,sessind1});
                        medianCorrMap(mapind,sessind1) = nanmedian(CorrMap{mapind,sessind1});
                        meanCorrMap(mapind,sessind1) = nanmean(CorrMap{mapind,sessind1});
                        [minCorrMap(mapind,sessind1), minrotCorrMap(mapind,sessind1)] = nanmin(CorrMap{mapind,sessind1});
                    end
                end
            end
        end                      
    end
end

%Cell count per session
out.Session1CellCount = length(find(Singlemap(:,1)));
out.Session2CellCount = length(find(Singlemap(:,2)));
out.Session3CellCount = length(find(Singlemap(:,3)));
out.Session4CellCount = length(find(Singlemap(:,4)));
%Map index of cells which passed passed EBC criteria
out.EBCind1 = ind1;
out.EBCind2 = ind2;
out.EBCind3 = ind3;
out.EBCind4 = ind4;
out.RMcorrMap = CorrMap;
out.RMcorrMapMax = maxCorrMap;
out.RMcorrMapMedian = medianCorrMap;
out.RMcorrMapMean = meanCorrMap;
out.RMcorrMapMin = minCorrMap;
out.RMcorrMapMaxRot = maxrotCorrMap;
out.RMcorrMapMinRot = minrotCorrMap;
end
