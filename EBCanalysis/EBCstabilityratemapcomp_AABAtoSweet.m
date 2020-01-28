%% Compares EBCs Ratemaps B context EBCs across A context EBCs
function [out] = EBCstabilityratemapcomp_AABAtoSweet(Singlemap,fullpass1,fullpass2,fullpass3, fullpass4, ebc1, ebc2, ebc3, ebc4, ebc5, ebc6, ebc7, ebc8, ebc9, ebc10)

[~, ind1, ~] = intersect(Singlemap(:,1),fullpass1);
[~, ind2, ~] = intersect(Singlemap(:,2),fullpass2);
[~, ind3, ~] = intersect(Singlemap(:,3),fullpass3);
[~, ind4, ~] = intersect(Singlemap(:,4),fullpass4);

EBCreg12 = Singlemap(ind1(find(Singlemap(ind1,2)>0)),2);
EBCreg21 = Singlemap(ind2(find(Singlemap(ind2,1)>0)),1);
EBCreg13 = Singlemap(ind1(find(Singlemap(ind1,3)>0)),3);
EBCreg31 = Singlemap(ind3(find(Singlemap(ind3,1)>0)),1);
EBCreg14 = Singlemap(ind1(find(Singlemap(ind1,4)>0)),4);
EBCreg41 = Singlemap(ind4(find(Singlemap(ind4,1)>0)),1);
EBCreg15 = Singlemap(ind1(find(Singlemap(ind1,5)>0)),5);
EBCreg16 = Singlemap(ind1(find(Singlemap(ind1,6)>0)),6);
EBCreg17 = Singlemap(ind1(find(Singlemap(ind1,7)>0)),7);
EBCreg18 = Singlemap(ind1(find(Singlemap(ind1,8)>0)),8);
EBCreg19 = Singlemap(ind1(find(Singlemap(ind1,9)>0)),9);
EBCreg110 = Singlemap(ind1(find(Singlemap(ind1,10)>0)),10);

EBCreg23 = Singlemap(ind2(find(Singlemap(ind2,3)>0)),3);
EBCreg24 = Singlemap(ind2(find(Singlemap(ind2,4)>0)),4);
EBCreg25 = Singlemap(ind2(find(Singlemap(ind2,5)>0)),5);
EBCreg26 = Singlemap(ind2(find(Singlemap(ind2,6)>0)),6);
EBCreg27 = Singlemap(ind2(find(Singlemap(ind2,7)>0)),7);
EBCreg28 = Singlemap(ind2(find(Singlemap(ind2,8)>0)),8);
EBCreg29 = Singlemap(ind2(find(Singlemap(ind2,9)>0)),9);
EBCreg210 = Singlemap(ind2(find(Singlemap(ind2,10)>0)),10);

EBCreg32 = Singlemap(ind3(find(Singlemap(ind3,2)>0)),2);
EBCreg34 = Singlemap(ind3(find(Singlemap(ind3,4)>0)),4);
EBCreg35 = Singlemap(ind3(find(Singlemap(ind3,5)>0)),5);
EBCreg36 = Singlemap(ind3(find(Singlemap(ind3,6)>0)),6);
EBCreg37 = Singlemap(ind3(find(Singlemap(ind3,7)>0)),7);
EBCreg38 = Singlemap(ind3(find(Singlemap(ind3,8)>0)),8);
EBCreg39 = Singlemap(ind3(find(Singlemap(ind3,9)>0)),9);
EBCreg310 = Singlemap(ind3(find(Singlemap(ind3,10)>0)),10);

EBCreg42 = Singlemap(ind4(find(Singlemap(ind4,2)>0)),2);
EBCreg43 = Singlemap(ind4(find(Singlemap(ind4,3)>0)),3);
EBCreg43 = Singlemap(ind4(find(Singlemap(ind4,3)>0)),3);
EBCreg45 = Singlemap(ind4(find(Singlemap(ind4,5)>0)),5);
EBCreg46 = Singlemap(ind4(find(Singlemap(ind4,6)>0)),6);
EBCreg47 = Singlemap(ind4(find(Singlemap(ind4,7)>0)),7);
EBCreg48 = Singlemap(ind4(find(Singlemap(ind4,8)>0)),8);
EBCreg49 = Singlemap(ind4(find(Singlemap(ind4,9)>0)),9);
EBCreg410 = Singlemap(ind4(find(Singlemap(ind4,10)>0)),10);

EBCind12 = ind1(find(Singlemap(ind1,2)>0));
EBCind21 = ind2(find(Singlemap(ind2,1)>0));
EBCind13 = ind1(find(Singlemap(ind1,3)>0));
EBCind31 = ind3(find(Singlemap(ind3,1)>0));
EBCind14 = ind1(find(Singlemap(ind1,4)>0));
EBCind41 = ind4(find(Singlemap(ind4,1)>0));
EBCind15 = ind1(find(Singlemap(ind1,5)>0));
EBCind16 = ind1(find(Singlemap(ind1,6)>0));
EBCind17 = ind1(find(Singlemap(ind1,7)>0));
EBCind18 = ind1(find(Singlemap(ind1,8)>0));
EBCind19 = ind1(find(Singlemap(ind1,9)>0));
EBCind110= ind1(find(Singlemap(ind1,10)>0));

EBCind23 = ind2(find(Singlemap(ind2,3)>0));
EBCind24 = ind2(find(Singlemap(ind2,4)>0));
EBCind25 = ind2(find(Singlemap(ind2,5)>0));
EBCind26 = ind2(find(Singlemap(ind2,6)>0));
EBCind27 = ind2(find(Singlemap(ind2,7)>0));
EBCind28 = ind2(find(Singlemap(ind2,8)>0));
EBCind29 = ind2(find(Singlemap(ind2,9)>0));
EBCind210= ind2(find(Singlemap(ind2,10)>0));

EBCind32 = ind3(find(Singlemap(ind3,2)>0));
EBCind34 = ind3(find(Singlemap(ind3,4)>0));
EBCind35 = ind3(find(Singlemap(ind3,5)>0));
EBCind36 = ind3(find(Singlemap(ind3,6)>0));
EBCind37 = ind3(find(Singlemap(ind3,7)>0));
EBCind38 = ind3(find(Singlemap(ind3,8)>0));
EBCind39 = ind3(find(Singlemap(ind3,9)>0));
EBCind310= ind3(find(Singlemap(ind3,10)>0));

EBCind42 = ind4(find(Singlemap(ind4,2)>0));
EBCind43 = ind4(find(Singlemap(ind4,3)>0));
EBCind45 = ind4(find(Singlemap(ind4,5)>0));
EBCind46 = ind4(find(Singlemap(ind4,6)>0));
EBCind47 = ind4(find(Singlemap(ind4,7)>0));
EBCind48 = ind4(find(Singlemap(ind4,8)>0));
EBCind49 = ind4(find(Singlemap(ind4,9)>0));
EBCind410= ind4(find(Singlemap(ind4,10)>0));

pass12 = intersect(ind1,ind2);
pass13 = intersect(ind1,ind3);
pass14 = intersect(ind1,ind4);
pass23 = intersect(ind2,ind3);
pass24 = intersect(ind2,ind4);
pass34 = intersect(ind3,ind4);

bmap = Singlemap;
bmap(bmap>0) = 1;

reg12 = find((bmap(:,1)+ bmap(:,2)) == 2);
reg13 = find((bmap(:,1)+ bmap(:,3)) == 2);
reg14 = find((bmap(:,1)+ bmap(:,4)) == 2);
reg15 = find((bmap(:,1)+ bmap(:,5)) == 2);
reg16 = find((bmap(:,1)+ bmap(:,6)) == 2);
reg17 = find((bmap(:,1)+ bmap(:,7)) == 2);
reg18 = find((bmap(:,1)+ bmap(:,8)) == 2);
reg19 = find((bmap(:,1)+ bmap(:,9)) == 2);
reg110= find((bmap(:,1)+ bmap(:,10))== 2);

reg23 = find((bmap(:,2)+ bmap(:,3)) == 2);
reg24 = find((bmap(:,2)+ bmap(:,4)) == 2);
reg25 = find((bmap(:,2)+ bmap(:,5)) == 2);
reg26 = find((bmap(:,2)+ bmap(:,6)) == 2);
reg27 = find((bmap(:,2)+ bmap(:,7)) == 2);
reg28 = find((bmap(:,2)+ bmap(:,8)) == 2);
reg29 = find((bmap(:,2)+ bmap(:,9)) == 2);
reg210=find((bmap(:,2)+ bmap(:,10)) == 2);

reg34 = find((bmap(:,4)+ bmap(:,3)) == 2);
reg35 = find((bmap(:,5)+ bmap(:,3)) == 2);
reg36 = find((bmap(:,6)+ bmap(:,3)) == 2);
reg37 = find((bmap(:,7)+ bmap(:,3)) == 2);
reg38 = find((bmap(:,8)+ bmap(:,3)) == 2);
reg39 = find((bmap(:,9)+ bmap(:,3)) == 2);
reg310= find((bmap(:,10)+bmap(:,3)) == 2);

reg45 = find((bmap(:,4)+ bmap(:,5)) == 2);
reg46 = find((bmap(:,4)+ bmap(:,6)) == 2);
reg47 = find((bmap(:,4)+ bmap(:,7)) == 2);
reg48 = find((bmap(:,4)+ bmap(:,8)) == 2);
reg49 = find((bmap(:,4)+ bmap(:,9)) == 2);
reg410= find((bmap(:,4)+bmap(:,10)) == 2);

reg56 = find((bmap(:,5)+ bmap(:,6)) == 2);
reg57 = find((bmap(:,5)+ bmap(:,7)) == 2);
reg58 = find((bmap(:,5)+ bmap(:,8)) == 2);
reg59 = find((bmap(:,5)+ bmap(:,9)) == 2);
reg510= find((bmap(:,5)+bmap(:,10)) == 2);

reg67 = find((bmap(:,6)+ bmap(:,7)) == 2);
reg68 = find((bmap(:,6)+ bmap(:,8)) == 2);
reg69 = find((bmap(:,6)+ bmap(:,9)) == 2);
reg610= find((bmap(:,6)+bmap(:,10)) == 2);

reg78 = find((bmap(:,7)+ bmap(:,8)) == 2);
reg79 = find((bmap(:,7)+ bmap(:,9)) == 2);
reg710= find((bmap(:,7)+bmap(:,10)) == 2);

reg89 = find((bmap(:,8)+ bmap(:,9)) == 2);
reg810= find((bmap(:,8)+bmap(:,10)) == 2);

reg910= find((bmap(:,9)+bmap(:,10)) == 2);

CorrMap = cell(size(Singlemap,1),size(Singlemap,2),10);
maxCorrMap = nan(size(Singlemap,1),size(Singlemap,2),10);
medianCorrMap = nan(size(Singlemap,1),size(Singlemap,2),10);
meanCorrMap = nan(size(Singlemap,1),size(Singlemap,2),10);
minCorrMap = nan(size(Singlemap,1),size(Singlemap,2),10);
maxrotCorrMap = nan(size(Singlemap,1),size(Singlemap,2),10);
minrotCorrMap = nan(size(Singlemap,1),size(Singlemap,2),10);

for mapind = 1 : length(Singlemap(:,1))
    if sum(bmap(mapind,:)) > 1
        for sessind1 = 1 : 9
            if Singlemap(mapind,sessind1) > 0
                if sessind1 == 1
                    if Singlemap(mapind,2) > 0
                        CorrMap{mapind,sessind1,2} = rmRotcorr_singlecell(ebc1.rm(:,:,Singlemap(mapind,1)),ebc2.rm(:,:,Singlemap(mapind,2)));
                        [maxCorrMap(mapind,sessind1,2), maxrotCorrMap(mapind,sessind1,2)] = nanmax(CorrMap{mapind,sessind1,2});
                        medianCorrMap(mapind,sessind1,2) = nanmedian(CorrMap{mapind,sessind1,2});
                        meanCorrMap(mapind,sessind1,2) = nanmean(CorrMap{mapind,sessind1,2});
                        [minCorrMap(mapind,sessind1,2), minrotCorrMap(mapind,sessind1,2)] = nanmin(CorrMap{mapind,sessind1,2});
                    end
                    if Singlemap(mapind,3) > 0
                        CorrMap{mapind,sessind1,3} = rmRotcorr_singlecell(ebc1.rm(:,:,Singlemap(mapind,1)),ebc3.rm(1:25,:,Singlemap(mapind,3)));
                        [maxCorrMap(mapind,sessind1,3), maxrotCorrMap(mapind,sessind1,3)] = nanmax(CorrMap{mapind,sessind1,3});
                        medianCorrMap(mapind,sessind1,3) = nanmedian(CorrMap{mapind,sessind1,3});
                        meanCorrMap(mapind,sessind1,3) = nanmean(CorrMap{mapind,sessind1,3});
                        [minCorrMap(mapind,sessind1,3), minrotCorrMap(mapind,sessind1,3)] = nanmin(CorrMap{mapind,sessind1,3});
                    end
                    if Singlemap(mapind,4) > 0
                        CorrMap{mapind,sessind1,4} = rmRotcorr_singlecell(ebc1.rm(:,:,Singlemap(mapind,1)),ebc4.rm(:,:,Singlemap(mapind,4)));
                        [maxCorrMap(mapind,sessind1,4), maxrotCorrMap(mapind,sessind1,4)] = nanmax(CorrMap{mapind,sessind1,4});
                        medianCorrMap(mapind,sessind1,4) = nanmedian(CorrMap{mapind,sessind1,4});
                        meanCorrMap(mapind,sessind1,4) = nanmean(CorrMap{mapind,sessind1,4});
                        [minCorrMap(mapind,sessind1,4), minrotCorrMap(mapind,sessind1,4)] = nanmin(CorrMap{mapind,sessind1,4});
                    end
                    if Singlemap(mapind,5) > 0
                        CorrMap{mapind,sessind1,5} = rmRotcorr_singlecell(ebc1.rm(:,:,Singlemap(mapind,1)),ebc5.rm(:,:,Singlemap(mapind,5)));
                        [maxCorrMap(mapind,sessind1,5), maxrotCorrMap(mapind,sessind1,5)] = nanmax(CorrMap{mapind,sessind1,5});
                        medianCorrMap(mapind,sessind1,5) = nanmedian(CorrMap{mapind,sessind1,5});
                        meanCorrMap(mapind,sessind1,5) = nanmean(CorrMap{mapind,sessind1,5});
                        [minCorrMap(mapind,sessind1,5), minrotCorrMap(mapind,sessind1,5)] = nanmin(CorrMap{mapind,sessind1,5});
                    end
                    if Singlemap(mapind,6) > 0
                        CorrMap{mapind,sessind1,6} = rmRotcorr_singlecell(ebc1.rm(:,:,Singlemap(mapind,1)),ebc6.rm(:,:,Singlemap(mapind,6)));
                        [maxCorrMap(mapind,sessind1,6), maxrotCorrMap(mapind,sessind1,6)] = nanmax(CorrMap{mapind,sessind1,6});
                        medianCorrMap(mapind,sessind1,6) = nanmedian(CorrMap{mapind,sessind1,6});
                        meanCorrMap(mapind,sessind1,6) = nanmean(CorrMap{mapind,sessind1,6});
                        [minCorrMap(mapind,sessind1,6), minrotCorrMap(mapind,sessind1,6)] = nanmin(CorrMap{mapind,sessind1,6});
                    end
                    if Singlemap(mapind,7) > 0
                        CorrMap{mapind,sessind1,7} = rmRotcorr_singlecell(ebc1.rm(:,:,Singlemap(mapind,1)),ebc7.rm(:,:,Singlemap(mapind,7)));
                        [maxCorrMap(mapind,sessind1,7), maxrotCorrMap(mapind,sessind1,7)] = nanmax(CorrMap{mapind,sessind1,7});
                        medianCorrMap(mapind,sessind1,7) = nanmedian(CorrMap{mapind,sessind1,7});
                        meanCorrMap(mapind,sessind1,7) = nanmean(CorrMap{mapind,sessind1,7});
                        [minCorrMap(mapind,sessind1,7), minrotCorrMap(mapind,sessind1,7)] = nanmin(CorrMap{mapind,sessind1,7});
                    end
                    if Singlemap(mapind,8) > 0
                        CorrMap{mapind,sessind1,8} = rmRotcorr_singlecell(ebc1.rm(:,:,Singlemap(mapind,1)),ebc8.rm(:,:,Singlemap(mapind,8)));
                        [maxCorrMap(mapind,sessind1,8), maxrotCorrMap(mapind,sessind1,8)] = nanmax(CorrMap{mapind,sessind1,8});
                        medianCorrMap(mapind,sessind1,8) = nanmedian(CorrMap{mapind,sessind1,8});
                        meanCorrMap(mapind,sessind1,8) = nanmean(CorrMap{mapind,sessind1,8});
                        [minCorrMap(mapind,sessind1,8), minrotCorrMap(mapind,sessind1,8)] = nanmin(CorrMap{mapind,sessind1,8});
                    end
                    if Singlemap(mapind,9) > 0
                        CorrMap{mapind,sessind1,9} = rmRotcorr_singlecell(ebc1.rm(:,:,Singlemap(mapind,1)),ebc9.rm(:,:,Singlemap(mapind,9)));
                        [maxCorrMap(mapind,sessind1,9), maxrotCorrMap(mapind,sessind1,9)] = nanmax(CorrMap{mapind,sessind1,9});
                        medianCorrMap(mapind,sessind1,9) = nanmedian(CorrMap{mapind,sessind1,9});
                        meanCorrMap(mapind,sessind1,9) = nanmean(CorrMap{mapind,sessind1,9});
                        [minCorrMap(mapind,sessind1,9), minrotCorrMap(mapind,sessind1,9)] = nanmin(CorrMap{mapind,sessind1,9});
                    end
                    if Singlemap(mapind,10) > 0
                        CorrMap{mapind,sessind1,10}= rmRotcorr_singlecell(ebc1.rm(:,:,Singlemap(mapind,1)),ebc10.rm(:,:,Singlemap(mapind,10)));
                        [maxCorrMap(mapind,sessind1,10), maxrotCorrMap(mapind,sessind1,10)] = nanmax(CorrMap{mapind,sessind1,10});
                        medianCorrMap(mapind,sessind1,10) = nanmedian(CorrMap{mapind,sessind1,10});
                        meanCorrMap(mapind,sessind1,10) = nanmean(CorrMap{mapind,sessind1,10});
                        [minCorrMap(mapind,sessind1,10), minrotCorrMap(mapind,sessind1,10)] = nanmin(CorrMap{mapind,sessind1,10});
                    end                    
                elseif sessind1 == 2
                    if Singlemap(mapind,3) > 0
                        CorrMap{mapind,sessind1,3} = rmRotcorr_singlecell(ebc2.rm(:,:,Singlemap(mapind,2)),ebc3.rm(1:25,:,Singlemap(mapind,3)));
                        [maxCorrMap(mapind,sessind1,3), maxrotCorrMap(mapind,sessind1,3)] = nanmax(CorrMap{mapind,sessind1,3});
                        medianCorrMap(mapind,sessind1,3) = nanmedian(CorrMap{mapind,sessind1,3});
                        meanCorrMap(mapind,sessind1,3) = nanmean(CorrMap{mapind,sessind1,3});
                        [minCorrMap(mapind,sessind1,3), minrotCorrMap(mapind,sessind1,3)] = nanmin(CorrMap{mapind,sessind1,3});
                    end
                    if Singlemap(mapind,4) > 0
                        CorrMap{mapind,sessind1,4} = rmRotcorr_singlecell(ebc2.rm(:,:,Singlemap(mapind,2)),ebc4.rm(:,:,Singlemap(mapind,4)));
                        [maxCorrMap(mapind,sessind1,4), maxrotCorrMap(mapind,sessind1,4)] = nanmax(CorrMap{mapind,sessind1,4});
                        medianCorrMap(mapind,sessind1,4) = nanmedian(CorrMap{mapind,sessind1,4});
                        meanCorrMap(mapind,sessind1,4) = nanmean(CorrMap{mapind,sessind1,4});
                        [minCorrMap(mapind,sessind1,4), minrotCorrMap(mapind,sessind1,4)] = nanmin(CorrMap{mapind,sessind1,4});
                    end
                    if Singlemap(mapind,5) > 0
                        CorrMap{mapind,sessind1,5} = rmRotcorr_singlecell(ebc2.rm(:,:,Singlemap(mapind,2)),ebc5.rm(:,:,Singlemap(mapind,5)));
                        [maxCorrMap(mapind,sessind1,5), maxrotCorrMap(mapind,sessind1,5)] = nanmax(CorrMap{mapind,sessind1,5});
                        medianCorrMap(mapind,sessind1,5) = nanmedian(CorrMap{mapind,sessind1,5});
                        meanCorrMap(mapind,sessind1,5) = nanmean(CorrMap{mapind,sessind1,5});
                        [minCorrMap(mapind,sessind1,5), minrotCorrMap(mapind,sessind1,5)] = nanmin(CorrMap{mapind,sessind1,5});
                    end
                    if Singlemap(mapind,6) > 0
                        CorrMap{mapind,sessind1,6} = rmRotcorr_singlecell(ebc2.rm(:,:,Singlemap(mapind,2)),ebc6.rm(:,:,Singlemap(mapind,6)));
                        [maxCorrMap(mapind,sessind1,6), maxrotCorrMap(mapind,sessind1,6)] = nanmax(CorrMap{mapind,sessind1,6});
                        medianCorrMap(mapind,sessind1,6) = nanmedian(CorrMap{mapind,sessind1,6});
                        meanCorrMap(mapind,sessind1,6) = nanmean(CorrMap{mapind,sessind1,6});
                        [minCorrMap(mapind,sessind1,6), minrotCorrMap(mapind,sessind1,6)] = nanmin(CorrMap{mapind,sessind1,6});
                    end
                    if Singlemap(mapind,7) > 0
                        CorrMap{mapind,sessind1,7} = rmRotcorr_singlecell(ebc2.rm(:,:,Singlemap(mapind,2)),ebc7.rm(:,:,Singlemap(mapind,7)));
                        [maxCorrMap(mapind,sessind1,7), maxrotCorrMap(mapind,sessind1,7)] = nanmax(CorrMap{mapind,sessind1,7});
                        medianCorrMap(mapind,sessind1,7) = nanmedian(CorrMap{mapind,sessind1,7});
                        meanCorrMap(mapind,sessind1,7) = nanmean(CorrMap{mapind,sessind1,7});
                        [minCorrMap(mapind,sessind1,7), minrotCorrMap(mapind,sessind1,7)] = nanmin(CorrMap{mapind,sessind1,7});
                    end
                    if Singlemap(mapind,8) > 0
                        CorrMap{mapind,sessind1,8} = rmRotcorr_singlecell(ebc2.rm(:,:,Singlemap(mapind,2)),ebc8.rm(:,:,Singlemap(mapind,8)));
                        [maxCorrMap(mapind,sessind1,8), maxrotCorrMap(mapind,sessind1,8)] = nanmax(CorrMap{mapind,sessind1,8});
                        medianCorrMap(mapind,sessind1,8) = nanmedian(CorrMap{mapind,sessind1,8});
                        meanCorrMap(mapind,sessind1,8) = nanmean(CorrMap{mapind,sessind1,8});
                        [minCorrMap(mapind,sessind1,8), minrotCorrMap(mapind,sessind1,8)] = nanmin(CorrMap{mapind,sessind1,8});
                    end
                    if Singlemap(mapind,9) > 0
                        CorrMap{mapind,sessind1,9} = rmRotcorr_singlecell(ebc2.rm(:,:,Singlemap(mapind,2)),ebc9.rm(:,:,Singlemap(mapind,9)));
                        [maxCorrMap(mapind,sessind1,9), maxrotCorrMap(mapind,sessind1,9)] = nanmax(CorrMap{mapind,sessind1,9});
                        medianCorrMap(mapind,sessind1,9) = nanmedian(CorrMap{mapind,sessind1,9});
                        meanCorrMap(mapind,sessind1,9) = nanmean(CorrMap{mapind,sessind1,9});
                        [minCorrMap(mapind,sessind1,9), minrotCorrMap(mapind,sessind1,9)] = nanmin(CorrMap{mapind,sessind1,9});
                    end
                    if Singlemap(mapind,10) > 0
                        CorrMap{mapind,sessind1,10}= rmRotcorr_singlecell(ebc2.rm(:,:,Singlemap(mapind,2)),ebc10.rm(:,:,Singlemap(mapind,10)));
                        [maxCorrMap(mapind,sessind1,10), maxrotCorrMap(mapind,sessind1,10)] = nanmax(CorrMap{mapind,sessind1,10});
                        medianCorrMap(mapind,sessind1,10) = nanmedian(CorrMap{mapind,sessind1,10});
                        meanCorrMap(mapind,sessind1,10) = nanmean(CorrMap{mapind,sessind1,10});
                        [minCorrMap(mapind,sessind1,10), minrotCorrMap(mapind,sessind1,10)] = nanmin(CorrMap{mapind,sessind1,10});
                    end 
                elseif sessind1 == 3
                    if Singlemap(mapind,4) > 0
                        CorrMap{mapind,sessind1,4} = rmRotcorr_singlecell(ebc3.rm(1:25,:,Singlemap(mapind,3)),ebc4.rm(:,:,Singlemap(mapind,4)));
                        [maxCorrMap(mapind,sessind1,4), maxrotCorrMap(mapind,sessind1,4)] = nanmax(CorrMap{mapind,sessind1,4});
                        medianCorrMap(mapind,sessind1,4) = nanmedian(CorrMap{mapind,sessind1,4});
                        meanCorrMap(mapind,sessind1,4) = nanmean(CorrMap{mapind,sessind1,4});
                        [minCorrMap(mapind,sessind1,4), minrotCorrMap(mapind,sessind1,4)] = nanmin(CorrMap{mapind,sessind1,4});
                    end
                    if Singlemap(mapind,5) > 0
                        CorrMap{mapind,sessind1,5} = rmRotcorr_singlecell(ebc3.rm(1:25,:,Singlemap(mapind,3)),ebc5.rm(:,:,Singlemap(mapind,5)));
                        [maxCorrMap(mapind,sessind1,5), maxrotCorrMap(mapind,sessind1,5)] = nanmax(CorrMap{mapind,sessind1,5});
                        medianCorrMap(mapind,sessind1,5) = nanmedian(CorrMap{mapind,sessind1,5});
                        meanCorrMap(mapind,sessind1,5) = nanmean(CorrMap{mapind,sessind1,5});
                        [minCorrMap(mapind,sessind1,5), minrotCorrMap(mapind,sessind1,5)] = nanmin(CorrMap{mapind,sessind1,5});
                    end
                    if Singlemap(mapind,6) > 0
                        CorrMap{mapind,sessind1,6} = rmRotcorr_singlecell(ebc3.rm(1:25,:,Singlemap(mapind,3)),ebc6.rm(:,:,Singlemap(mapind,6)));
                        [maxCorrMap(mapind,sessind1,6), maxrotCorrMap(mapind,sessind1,6)] = nanmax(CorrMap{mapind,sessind1,6});
                        medianCorrMap(mapind,sessind1,6) = nanmedian(CorrMap{mapind,sessind1,6});
                        meanCorrMap(mapind,sessind1,6) = nanmean(CorrMap{mapind,sessind1,6});
                        [minCorrMap(mapind,sessind1,6), minrotCorrMap(mapind,sessind1,6)] = nanmin(CorrMap{mapind,sessind1,6});
                    end
                    if Singlemap(mapind,7) > 0
                        CorrMap{mapind,sessind1,7} = rmRotcorr_singlecell(ebc3.rm(1:25,:,Singlemap(mapind,3)),ebc7.rm(:,:,Singlemap(mapind,7)));
                        [maxCorrMap(mapind,sessind1,7), maxrotCorrMap(mapind,sessind1,7)] = nanmax(CorrMap{mapind,sessind1,7});
                        medianCorrMap(mapind,sessind1,7) = nanmedian(CorrMap{mapind,sessind1,7});
                        meanCorrMap(mapind,sessind1,7) = nanmean(CorrMap{mapind,sessind1,7});
                        [minCorrMap(mapind,sessind1,7), minrotCorrMap(mapind,sessind1,7)] = nanmin(CorrMap{mapind,sessind1,7});
                    end
                    if Singlemap(mapind,8) > 0
                        CorrMap{mapind,sessind1,8} = rmRotcorr_singlecell(ebc3.rm(1:25,:,Singlemap(mapind,3)),ebc8.rm(:,:,Singlemap(mapind,8)));
                        [maxCorrMap(mapind,sessind1,8), maxrotCorrMap(mapind,sessind1,8)] = nanmax(CorrMap{mapind,sessind1,8});
                        medianCorrMap(mapind,sessind1,8) = nanmedian(CorrMap{mapind,sessind1,8});
                        meanCorrMap(mapind,sessind1,8) = nanmean(CorrMap{mapind,sessind1,8});
                        [minCorrMap(mapind,sessind1,8), minrotCorrMap(mapind,sessind1,8)] = nanmin(CorrMap{mapind,sessind1,8});
                    end
                    if Singlemap(mapind,9) > 0
                        CorrMap{mapind,sessind1,9} = rmRotcorr_singlecell(ebc3.rm(1:25,:,Singlemap(mapind,3)),ebc9.rm(:,:,Singlemap(mapind,9)));
                        [maxCorrMap(mapind,sessind1,9), maxrotCorrMap(mapind,sessind1,9)] = nanmax(CorrMap{mapind,sessind1,9});
                        medianCorrMap(mapind,sessind1,9) = nanmedian(CorrMap{mapind,sessind1,9});
                        meanCorrMap(mapind,sessind1,9) = nanmean(CorrMap{mapind,sessind1,9});
                        [minCorrMap(mapind,sessind1,9), minrotCorrMap(mapind,sessind1,9)] = nanmin(CorrMap{mapind,sessind1,9});
                    end
                    if Singlemap(mapind,10) > 0
                        CorrMap{mapind,sessind1,10}= rmRotcorr_singlecell(ebc3.rm(1:25,:,Singlemap(mapind,3)),ebc10.rm(:,:,Singlemap(mapind,10)));
                        [maxCorrMap(mapind,sessind1,10), maxrotCorrMap(mapind,sessind1,10)] = nanmax(CorrMap{mapind,sessind1,10});
                        medianCorrMap(mapind,sessind1,10) = nanmedian(CorrMap{mapind,sessind1,10});
                        meanCorrMap(mapind,sessind1,10) = nanmean(CorrMap{mapind,sessind1,10});
                        [minCorrMap(mapind,sessind1,10), minrotCorrMap(mapind,sessind1,10)] = nanmin(CorrMap{mapind,sessind1,10});
                    end 
                elseif sessind1 == 4
                    if Singlemap(mapind,5) > 0
                        CorrMap{mapind,sessind1,5} = rmRotcorr_singlecell(ebc4.rm(:,:,Singlemap(mapind,4)),ebc5.rm(:,:,Singlemap(mapind,5)));
                        [maxCorrMap(mapind,sessind1,5), maxrotCorrMap(mapind,sessind1,5)] = nanmax(CorrMap{mapind,sessind1,5});
                        medianCorrMap(mapind,sessind1,5) = nanmedian(CorrMap{mapind,sessind1,5});
                        meanCorrMap(mapind,sessind1,5) = nanmean(CorrMap{mapind,sessind1,5});
                        [minCorrMap(mapind,sessind1,5), minrotCorrMap(mapind,sessind1,5)] = nanmin(CorrMap{mapind,sessind1,5});
                    end
                    if Singlemap(mapind,6) > 0
                        CorrMap{mapind,sessind1,6} = rmRotcorr_singlecell(ebc4.rm(:,:,Singlemap(mapind,4)),ebc6.rm(:,:,Singlemap(mapind,6)));
                        [maxCorrMap(mapind,sessind1,6), maxrotCorrMap(mapind,sessind1,6)] = nanmax(CorrMap{mapind,sessind1,6});
                        medianCorrMap(mapind,sessind1,6) = nanmedian(CorrMap{mapind,sessind1,6});
                        meanCorrMap(mapind,sessind1,6) = nanmean(CorrMap{mapind,sessind1,6});
                        [minCorrMap(mapind,sessind1,6), minrotCorrMap(mapind,sessind1,6)] = nanmin(CorrMap{mapind,sessind1,6});
                    end
                    if Singlemap(mapind,7) > 0
                        CorrMap{mapind,sessind1,7} = rmRotcorr_singlecell(ebc4.rm(:,:,Singlemap(mapind,4)),ebc7.rm(:,:,Singlemap(mapind,7)));
                        [maxCorrMap(mapind,sessind1,7), maxrotCorrMap(mapind,sessind1,7)] = nanmax(CorrMap{mapind,sessind1,7});
                        medianCorrMap(mapind,sessind1,7) = nanmedian(CorrMap{mapind,sessind1,7});
                        meanCorrMap(mapind,sessind1,7) = nanmean(CorrMap{mapind,sessind1,7});
                        [minCorrMap(mapind,sessind1,7), minrotCorrMap(mapind,sessind1,7)] = nanmin(CorrMap{mapind,sessind1,7});
                    end
                    if Singlemap(mapind,8) > 0
                        CorrMap{mapind,sessind1,8} = rmRotcorr_singlecell(ebc4.rm(:,:,Singlemap(mapind,4)),ebc8.rm(:,:,Singlemap(mapind,8)));
                        [maxCorrMap(mapind,sessind1,8), maxrotCorrMap(mapind,sessind1,8)] = nanmax(CorrMap{mapind,sessind1,8});
                        medianCorrMap(mapind,sessind1,8) = nanmedian(CorrMap{mapind,sessind1,8});
                        meanCorrMap(mapind,sessind1,8) = nanmean(CorrMap{mapind,sessind1,8});
                        [minCorrMap(mapind,sessind1,8), minrotCorrMap(mapind,sessind1,8)] = nanmin(CorrMap{mapind,sessind1,8});
                    end
                    if Singlemap(mapind,9) > 0
                        CorrMap{mapind,sessind1,9} = rmRotcorr_singlecell(ebc4.rm(:,:,Singlemap(mapind,4)),ebc9.rm(:,:,Singlemap(mapind,9)));
                        [maxCorrMap(mapind,sessind1,9), maxrotCorrMap(mapind,sessind1,9)] = nanmax(CorrMap{mapind,sessind1,9});
                        medianCorrMap(mapind,sessind1,9) = nanmedian(CorrMap{mapind,sessind1,9});
                        meanCorrMap(mapind,sessind1,9) = nanmean(CorrMap{mapind,sessind1,9});
                        [minCorrMap(mapind,sessind1,9), minrotCorrMap(mapind,sessind1,9)] = nanmin(CorrMap{mapind,sessind1,9});
                    end
                    if Singlemap(mapind,10) > 0
                        CorrMap{mapind,sessind1,10}= rmRotcorr_singlecell(ebc4.rm(:,:,Singlemap(mapind,4)),ebc10.rm(:,:,Singlemap(mapind,10)));
                        [maxCorrMap(mapind,sessind1,10), maxrotCorrMap(mapind,sessind1,10)] = nanmax(CorrMap{mapind,sessind1,10});
                        medianCorrMap(mapind,sessind1,10) = nanmedian(CorrMap{mapind,sessind1,10});
                        meanCorrMap(mapind,sessind1,10) = nanmean(CorrMap{mapind,sessind1,10});
                        [minCorrMap(mapind,sessind1,10), minrotCorrMap(mapind,sessind1,10)] = nanmin(CorrMap{mapind,sessind1,10});
                    end 
                elseif sessind1 == 5                    
                    if Singlemap(mapind,6) > 0
                        CorrMap{mapind,sessind1,6} = rmRotcorr_singlecell(ebc5.rm(:,:,Singlemap(mapind,5)),ebc6.rm(:,:,Singlemap(mapind,6)));
                        [maxCorrMap(mapind,sessind1,6), maxrotCorrMap(mapind,sessind1,6)] = nanmax(CorrMap{mapind,sessind1,6});
                        medianCorrMap(mapind,sessind1,6) = nanmedian(CorrMap{mapind,sessind1,6});
                        meanCorrMap(mapind,sessind1,6) = nanmean(CorrMap{mapind,sessind1,6});
                        [minCorrMap(mapind,sessind1,6), minrotCorrMap(mapind,sessind1,6)] = nanmin(CorrMap{mapind,sessind1,6});
                    end
                    if Singlemap(mapind,7) > 0
                        CorrMap{mapind,sessind1,7} = rmRotcorr_singlecell(ebc5.rm(:,:,Singlemap(mapind,5)),ebc7.rm(:,:,Singlemap(mapind,7)));
                        [maxCorrMap(mapind,sessind1,7), maxrotCorrMap(mapind,sessind1,7)] = nanmax(CorrMap{mapind,sessind1,7});
                        medianCorrMap(mapind,sessind1,7) = nanmedian(CorrMap{mapind,sessind1,7});
                        meanCorrMap(mapind,sessind1,7) = nanmean(CorrMap{mapind,sessind1,7});
                        [minCorrMap(mapind,sessind1,7), minrotCorrMap(mapind,sessind1,7)] = nanmin(CorrMap{mapind,sessind1,7});
                    end
                    if Singlemap(mapind,8) > 0
                        CorrMap{mapind,sessind1,8} = rmRotcorr_singlecell(ebc5.rm(:,:,Singlemap(mapind,5)),ebc8.rm(:,:,Singlemap(mapind,8)));
                        [maxCorrMap(mapind,sessind1,8), maxrotCorrMap(mapind,sessind1,8)] = nanmax(CorrMap{mapind,sessind1,8});
                        medianCorrMap(mapind,sessind1,8) = nanmedian(CorrMap{mapind,sessind1,8});
                        meanCorrMap(mapind,sessind1,8) = nanmean(CorrMap{mapind,sessind1,8});
                        [minCorrMap(mapind,sessind1,8), minrotCorrMap(mapind,sessind1,8)] = nanmin(CorrMap{mapind,sessind1,8});
                    end
                    if Singlemap(mapind,9) > 0
                        CorrMap{mapind,sessind1,9} = rmRotcorr_singlecell(ebc5.rm(:,:,Singlemap(mapind,5)),ebc9.rm(:,:,Singlemap(mapind,9)));
                        [maxCorrMap(mapind,sessind1,9), maxrotCorrMap(mapind,sessind1,9)] = nanmax(CorrMap{mapind,sessind1,9});
                        medianCorrMap(mapind,sessind1,9) = nanmedian(CorrMap{mapind,sessind1,9});
                        meanCorrMap(mapind,sessind1,9) = nanmean(CorrMap{mapind,sessind1,9});
                        [minCorrMap(mapind,sessind1,9), minrotCorrMap(mapind,sessind1,9)] = nanmin(CorrMap{mapind,sessind1,9});
                    end
                    if Singlemap(mapind,10) > 0
                        CorrMap{mapind,sessind1,10}= rmRotcorr_singlecell(ebc5.rm(:,:,Singlemap(mapind,5)),ebc10.rm(:,:,Singlemap(mapind,10)));
                        [maxCorrMap(mapind,sessind1,10), maxrotCorrMap(mapind,sessind1,10)] = nanmax(CorrMap{mapind,sessind1,10});
                        medianCorrMap(mapind,sessind1,10) = nanmedian(CorrMap{mapind,sessind1,10});
                        meanCorrMap(mapind,sessind1,10) = nanmean(CorrMap{mapind,sessind1,10});
                        [minCorrMap(mapind,sessind1,10), minrotCorrMap(mapind,sessind1,10)] = nanmin(CorrMap{mapind,sessind1,10});
                    end 
                elseif sessind1 == 6                                        
                    if Singlemap(mapind,7) > 0
                        CorrMap{mapind,sessind1,7} = rmRotcorr_singlecell(ebc6.rm(:,:,Singlemap(mapind,6)),ebc7.rm(:,:,Singlemap(mapind,7)));
                        [maxCorrMap(mapind,sessind1,7), maxrotCorrMap(mapind,sessind1,7)] = nanmax(CorrMap{mapind,sessind1,7});
                        medianCorrMap(mapind,sessind1,7) = nanmedian(CorrMap{mapind,sessind1,7});
                        meanCorrMap(mapind,sessind1,7) = nanmean(CorrMap{mapind,sessind1,7});
                        [minCorrMap(mapind,sessind1,7), minrotCorrMap(mapind,sessind1,7)] = nanmin(CorrMap{mapind,sessind1,7});
                    end
                    if Singlemap(mapind,8) > 0
                        CorrMap{mapind,sessind1,8} = rmRotcorr_singlecell(ebc6.rm(:,:,Singlemap(mapind,6)),ebc8.rm(:,:,Singlemap(mapind,8)));
                        [maxCorrMap(mapind,sessind1,8), maxrotCorrMap(mapind,sessind1,8)] = nanmax(CorrMap{mapind,sessind1,8});
                        medianCorrMap(mapind,sessind1,8) = nanmedian(CorrMap{mapind,sessind1,8});
                        meanCorrMap(mapind,sessind1,8) = nanmean(CorrMap{mapind,sessind1,8});
                        [minCorrMap(mapind,sessind1,8), minrotCorrMap(mapind,sessind1,8)] = nanmin(CorrMap{mapind,sessind1,8});
                    end
                    if Singlemap(mapind,9) > 0
                        CorrMap{mapind,sessind1,9} = rmRotcorr_singlecell(ebc6.rm(:,:,Singlemap(mapind,6)),ebc9.rm(:,:,Singlemap(mapind,9)));
                        [maxCorrMap(mapind,sessind1,9), maxrotCorrMap(mapind,sessind1,9)] = nanmax(CorrMap{mapind,sessind1,9});
                        medianCorrMap(mapind,sessind1,9) = nanmedian(CorrMap{mapind,sessind1,9});
                        meanCorrMap(mapind,sessind1,9) = nanmean(CorrMap{mapind,sessind1,9});
                        [minCorrMap(mapind,sessind1,9), minrotCorrMap(mapind,sessind1,9)] = nanmin(CorrMap{mapind,sessind1,9});
                    end
                    if Singlemap(mapind,10) > 0
                        CorrMap{mapind,sessind1,10}= rmRotcorr_singlecell(ebc6.rm(:,:,Singlemap(mapind,6)),ebc10.rm(:,:,Singlemap(mapind,10)));
                        [maxCorrMap(mapind,sessind1,10), maxrotCorrMap(mapind,sessind1,10)] = nanmax(CorrMap{mapind,sessind1,10});
                        medianCorrMap(mapind,sessind1,10) = nanmedian(CorrMap{mapind,sessind1,10});
                        meanCorrMap(mapind,sessind1,10) = nanmean(CorrMap{mapind,sessind1,10});
                        [minCorrMap(mapind,sessind1,10), minrotCorrMap(mapind,sessind1,10)] = nanmin(CorrMap{mapind,sessind1,10});
                    end 
                elseif sessind1 == 7                    
                    if Singlemap(mapind,8) > 0
                        CorrMap{mapind,sessind1,8} = rmRotcorr_singlecell(ebc7.rm(:,:,Singlemap(mapind,7)),ebc8.rm(:,:,Singlemap(mapind,8)));
                        [maxCorrMap(mapind,sessind1,8), maxrotCorrMap(mapind,sessind1,8)] = nanmax(CorrMap{mapind,sessind1,8});
                        medianCorrMap(mapind,sessind1,8) = nanmedian(CorrMap{mapind,sessind1,8});
                        meanCorrMap(mapind,sessind1,8) = nanmean(CorrMap{mapind,sessind1,8});
                        [minCorrMap(mapind,sessind1,8), minrotCorrMap(mapind,sessind1,8)] = nanmin(CorrMap{mapind,sessind1,8});
                    end
                    if Singlemap(mapind,9) > 0
                        CorrMap{mapind,sessind1,9} = rmRotcorr_singlecell(ebc7.rm(:,:,Singlemap(mapind,7)),ebc9.rm(:,:,Singlemap(mapind,9)));
                        [maxCorrMap(mapind,sessind1,9), maxrotCorrMap(mapind,sessind1,9)] = nanmax(CorrMap{mapind,sessind1,9});
                        medianCorrMap(mapind,sessind1,9) = nanmedian(CorrMap{mapind,sessind1,9});
                        meanCorrMap(mapind,sessind1,9) = nanmean(CorrMap{mapind,sessind1,9});
                        [minCorrMap(mapind,sessind1,9), minrotCorrMap(mapind,sessind1,9)] = nanmin(CorrMap{mapind,sessind1,9});
                    end
                    if Singlemap(mapind,10) > 0
                        CorrMap{mapind,sessind1,10}= rmRotcorr_singlecell(ebc7.rm(:,:,Singlemap(mapind,7)),ebc10.rm(:,:,Singlemap(mapind,10)));
                        [maxCorrMap(mapind,sessind1,10), maxrotCorrMap(mapind,sessind1,10)] = nanmax(CorrMap{mapind,sessind1,10});
                        medianCorrMap(mapind,sessind1,10) = nanmedian(CorrMap{mapind,sessind1,10});
                        meanCorrMap(mapind,sessind1,10) = nanmean(CorrMap{mapind,sessind1,10});
                        [minCorrMap(mapind,sessind1,10), minrotCorrMap(mapind,sessind1,10)] = nanmin(CorrMap{mapind,sessind1,10});
                    end 
                elseif sessind1 == 8
                    if Singlemap(mapind,9) > 0                        
                        CorrMap{mapind,sessind1,9} = rmRotcorr_singlecell(ebc8.rm(:,:,Singlemap(mapind,8)),ebc9.rm(:,:,Singlemap(mapind,9)));
                        [maxCorrMap(mapind,sessind1,9), maxrotCorrMap(mapind,sessind1,9)] = nanmax(CorrMap{mapind,sessind1,9});
                        medianCorrMap(mapind,sessind1,9) = nanmedian(CorrMap{mapind,sessind1,9});
                        meanCorrMap(mapind,sessind1,9) = nanmean(CorrMap{mapind,sessind1,9});
                        [minCorrMap(mapind,sessind1,9), minrotCorrMap(mapind,sessind1,9)] = nanmin(CorrMap{mapind,sessind1,9});
                    end
                    if Singlemap(mapind,10) > 0
                        CorrMap{mapind,sessind1,10}= rmRotcorr_singlecell(ebc8.rm(:,:,Singlemap(mapind,8)),ebc10.rm(:,:,Singlemap(mapind,10)));
                        [maxCorrMap(mapind,sessind1,10), maxrotCorrMap(mapind,sessind1,10)] = nanmax(CorrMap{mapind,sessind1,10});
                        medianCorrMap(mapind,sessind1,10) = nanmedian(CorrMap{mapind,sessind1,10});
                        meanCorrMap(mapind,sessind1,10) = nanmean(CorrMap{mapind,sessind1,10});
                        [minCorrMap(mapind,sessind1,10), minrotCorrMap(mapind,sessind1,10)] = nanmin(CorrMap{mapind,sessind1,10});
                    end 
                elseif sessind1 == 9                    
                    if Singlemap(mapind,10) > 0
                        CorrMap{mapind,sessind1,10}= rmRotcorr_singlecell(ebc9.rm(:,:,Singlemap(mapind,9)),ebc10.rm(:,:,Singlemap(mapind,10)));
                        [maxCorrMap(mapind,sessind1,10), maxrotCorrMap(mapind,sessind1,10)] = nanmax(CorrMap{mapind,sessind1,10});
                        medianCorrMap(mapind,sessind1,10) = nanmedian(CorrMap{mapind,sessind1,10});
                        meanCorrMap(mapind,sessind1,10) = nanmean(CorrMap{mapind,sessind1,10});
                        [minCorrMap(mapind,sessind1,10), minrotCorrMap(mapind,sessind1,10)] = nanmin(CorrMap{mapind,sessind1,10});
                    end 
                end
            end
        end                      
    end
end

%Cells which pass criteria in both sessions
out.EBCkeep12 = pass12;
out.EBCkeep13 = pass13;
out.EBCkeep14 = pass14;
out.EBCkeep23 = pass23;
out.EBCkeep24 = pass24;
out.EBCkeep43 = pass34;

%Cells registered across sessions
out.RegisteredMapIndex12 = reg12;
out.RegisteredMapIndex13 = reg13;
out.RegisteredMapIndex14 = reg14;
out.RegisteredMapIndex15 = reg15;
out.RegisteredMapIndex16 = reg16;
out.RegisteredMapIndex17 = reg17;
out.RegisteredMapIndex18 = reg18;
out.RegisteredMapIndex19 = reg19;
out.RegisteredMapIndex110 = reg110;

out.RegisteredMapIndex23 = reg23;
out.RegisteredMapIndex24 = reg24;
out.RegisteredMapIndex25 = reg25;
out.RegisteredMapIndex26 = reg26;
out.RegisteredMapIndex27 = reg27;
out.RegisteredMapIndex28 = reg28;
out.RegisteredMapIndex29 = reg29;
out.RegisteredMapIndex210 = reg210;

out.RegisteredMapIndex34 = reg34;
out.RegisteredMapIndex35 = reg35;
out.RegisteredMapIndex36 = reg36;
out.RegisteredMapIndex37 = reg37;
out.RegisteredMapIndex38 = reg38;
out.RegisteredMapIndex39 = reg39;
out.RegisteredMapIndex310 = reg310;

out.RegisteredMapIndex45 = reg45;
out.RegisteredMapIndex46 = reg46;
out.RegisteredMapIndex47 = reg47;
out.RegisteredMapIndex48 = reg48;
out.RegisteredMapIndex49 = reg49;
out.RegisteredMapIndex410 = reg410;

out.RegisteredMapIndex56 = reg56;
out.RegisteredMapIndex57 = reg57;
out.RegisteredMapIndex58 = reg58;
out.RegisteredMapIndex59 = reg59;
out.RegisteredMapIndex510 = reg510;

out.RegisteredMapIndex67 = reg67;
out.RegisteredMapIndex68 = reg68;
out.RegisteredMapIndex69 = reg69;
out.RegisteredMapIndex610 = reg610;

out.RegisteredMapIndex78 = reg78;
out.RegisteredMapIndex79 = reg79;
out.RegisteredMapIndex710 = reg710;

out.RegisteredMapIndex89 = reg89;
out.RegisteredMapIndex810 = reg810;

out.RegisteredMapIndex910 = reg910;

%Cell count per session
out.Session1CellCount = length(find(Singlemap(:,1)));
out.Session2CellCount = length(find(Singlemap(:,2)));
out.Session3CellCount = length(find(Singlemap(:,3)));
out.Session4CellCount = length(find(Singlemap(:,4)));
out.Session5CellCount = length(find(Singlemap(:,5)));
out.Session6CellCount = length(find(Singlemap(:,6)));
out.Session7CellCount = length(find(Singlemap(:,7)));
out.Session8CellCount = length(find(Singlemap(:,8)));
out.Session9CellCount = length(find(Singlemap(:,9)));
out.Session10CellCount = length(find(Singlemap(:,10)));

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

%{
%Cells which were registered across 
out.EBC1reg2 = EBCreg12;
out.EBC2reg1 = EBCreg21;
out.EBC1reg3 = EBCreg13;
out.EBC3reg1 = EBCreg31;
out.EBC1reg4 = EBCreg14;
out.EBC4reg1 = EBCreg41;
out.EBC1reg5 = EBCreg15;
out.EBC1reg6 = EBCreg16;
out.EBC1reg7 = EBCreg17;
out.EBC1reg8 = EBCreg18;
out.EBC1reg9 = EBCreg19;
out.EBC1reg10 = EBCreg110;

out.EBC2reg3 = EBCreg23;
out.EBC3reg2 = EBCreg32;
out.EBC2reg4 = EBCreg24;
out.EBC2reg5 = EBCreg25;
out.EBC2reg6 = EBCreg26;
out.EBC2reg7 = EBCreg27;
out.EBC2reg8 = EBCreg28;
out.EBC2reg9 = EBCreg29;
out.EBC2reg10 = EBCreg210;

out.EBC3reg4 = EBCreg34;
out.EBC3reg5 = EBCreg35;
out.EBC3reg6 = EBCreg36;
out.EBC3reg7 = EBCreg37;
out.EBC3reg8 = EBCreg38;
out.EBC3reg9 = EBCreg39;
out.EBC3reg10 = EBCreg310;

out.EBC4reg2 = EBCreg42;
out.EBC4reg3 = EBCreg43;
out.EBC4reg5 = EBCreg45;
out.EBC4reg6 = EBCreg46;
out.EBC4reg7 = EBCreg47;
out.EBC4reg8 = EBCreg48;
out.EBC4reg9 = EBCreg49;
out.EBC4reg10 = EBCreg410;

out.EBCmapind12 = EBCind12;
out.EBCmapind21 = EBCind21;
out.EBCmapind13 = EBCind13;
out.EBCmapind31 = EBCind31;
out.EBCmapind14 = EBCind14;
out.EBCmapind41 = EBCind41;
out.EBCmapind15 = EBCind15;
out.EBCmapind16 = EBCind16;
out.EBCmapind17 = EBCind17;
out.EBCmapind18 = EBCind18;
out.EBCmapind19 = EBCind19;
out.EBCmapind110 = EBCind110;
out.EBCmapind23 = EBCind23;
out.EBCmapind32 = EBCind32;
out.EBCmapind24 = EBCind24;
out.EBCmapind25 = EBCind25;
out.EBCmapind26 = EBCind26;
out.EBCmapind27 = EBCind27;
out.EBCmapind28 = EBCind28;
out.EBCmapind29 = EBCind29;
out.EBCmapind210 = EBCind210;

out.EBCmapind34 = EBCind34;
out.EBCmapind35 = EBCind35;
out.EBCmapind36 = EBCind36;
out.EBCmapind37 = EBCind37;
out.EBCmapind38 = EBCind38;
out.EBCmapind39 = EBCind39;
out.EBCmapind310 = EBCind310;

out.EBCmapind42 = EBCind42;
out.EBCmapind43 = EBCind43;
out.EBCmapind45 = EBCind45;
out.EBCmapind46 = EBCind46;
out.EBCmapind47 = EBCind47;
out.EBCmapind48 = EBCind48;
out.EBCmapind49 = EBCind49;
out.EBCmapind410 = EBCind410;
%}
end
%{
%EBC ratemap correlation
for i = 1 : length(reg12)
    out.rmcor12(i) = corr2(ebc1.rm(:,:,Singlemap(reg12(i),1)),ebc2.rm(:,:,Singlemap(reg12(i),2)));
end
for i = 1 : length(reg13)
    out.rmcor13(i) = corr2(ebc1.rm(:,:,Singlemap(reg13(i),1)),ebc3.rm(1:25,:,Singlemap(reg13(i),3)));
end
for i = 1 : length(reg14)
    out.rmcor14(i) = corr2(ebc1.rm(:,:,Singlemap(reg14(i),1)),ebc4.rm(:,:,Singlemap(reg14(i),4)));
end
for i = 1 : length(reg15)
    out.rmcor15(i) = corr2(ebc1.rm(:,:,Singlemap(reg15(i),1)),ebc5.rm(:,:,Singlemap(reg15(i),5)));
end
for i = 1 : length(reg16)
    out.rmcor16(i) = corr2(ebc1.rm(:,:,Singlemap(reg16(i),1)),ebc6.rm(:,:,Singlemap(reg16(i),6)));
end
for i = 1 : length(reg17)
    out.rmcor17(i) = corr2(ebc1.rm(:,:,Singlemap(reg17(i),1)),ebc7.rm(:,:,Singlemap(reg17(i),7)));
end
for i = 1 : length(reg18)
    out.rmcor18(i) = corr2(ebc1.rm(:,:,Singlemap(reg18(i),1)),ebc8.rm(:,:,Singlemap(reg18(i),8)));
end
for i = 1 : length(reg19)
    out.rmcor19(i) = corr2(ebc1.rm(:,:,Singlemap(reg19(i),1)),ebc9.rm(:,:,Singlemap(reg19(i),9)));
end
for i = 1 : length(reg110)
    out.rmcor110(i) = corr2(ebc1.rm(:,:,Singlemap(reg110(i),1)),ebc10.rm(:,:,Singlemap(reg110(i),10)));
end

for i = 1 : length(reg23)
    out.rmcor23(i) = corr2(ebc2.rm(:,:,Singlemap(reg23(i),2)),ebc3.rm(1:25,:,Singlemap(reg23(i),3)));
end
for i = 1 : length(reg24)
    out.rmcor24(i) = corr2(ebc2.rm(:,:,Singlemap(reg24(i),2)),ebc4.rm(:,:,Singlemap(reg24(i),4)));
end
for i = 1 : length(reg25)
    out.rmcor25(i) = corr2(ebc2.rm(:,:,Singlemap(reg25(i),2)),ebc5.rm(:,:,Singlemap(reg25(i),5)));
end
for i = 1 : length(reg26)
    out.rmcor26(i) = corr2(ebc2.rm(:,:,Singlemap(reg26(i),2)),ebc6.rm(:,:,Singlemap(reg26(i),6)));
end
for i = 1 : length(reg27)
    out.rmcor27(i) = corr2(ebc2.rm(:,:,Singlemap(reg27(i),2)),ebc7.rm(:,:,Singlemap(reg27(i),7)));
end
for i = 1 : length(reg28)
    out.rmcor28(i) = corr2(ebc2.rm(:,:,Singlemap(reg28(i),2)),ebc8.rm(:,:,Singlemap(reg28(i),8)));
end
for i = 1 : length(reg29)
    out.rmcor29(i) = corr2(ebc2.rm(:,:,Singlemap(reg29(i),2)),ebc9.rm(:,:,Singlemap(reg29(i),9)));
end
for i = 1 : length(reg210)
    out.rmcor210(i) = corr2(ebc2.rm(:,:,Singlemap(reg210(i),2)),ebc10.rm(:,:,Singlemap(reg210(i),10)));
end


for i = 1 : length(reg34)
    out.rmcor34(i) = corr2(ebc4.rm(:,:,Singlemap(reg34(i),4)),ebc3.rm(1:25,:,Singlemap(reg34(i),3)));
end
for i = 1 : length(reg35)
    out.rmcor35(i) = corr2(ebc5.rm(:,:,Singlemap(reg35(i),5)),ebc3.rm(1:25,:,Singlemap(reg35(i),3)));
end
for i = 1 : length(reg36)
    out.rmcor36(i) = corr2(ebc6.rm(:,:,Singlemap(reg36(i),6)),ebc3.rm(1:25,:,Singlemap(reg36(i),3)));
end
for i = 1 : length(reg37)
    out.rmcor37(i) = corr2(ebc7.rm(:,:,Singlemap(reg34(i),7)),ebc3.rm(1:25,:,Singlemap(reg37(i),3)));
end
for i = 1 : length(reg38)
    out.rmcor38(i) = corr2(ebc8.rm(:,:,Singlemap(reg38(i),8)),ebc3.rm(1:25,:,Singlemap(reg38(i),3)));
end
for i = 1 : length(reg39)
    out.rmcor39(i) = corr2(ebc9.rm(:,:,Singlemap(reg39(i),4)),ebc3.rm(1:25,:,Singlemap(reg39(i),3)));
end
for i = 1 : length(reg310)
    out.rmcor310(i) = corr2(ebc10.rm(:,:,Singlemap(reg310(i),10)),ebc3.rm(1:25,:,Singlemap(reg310(i),3)));
end

for i = 1 : length(reg45)
    out.rmcor45(i) = corr2(ebc4.rm(:,:,Singlemap(reg45(i),4)),ebc5.rm(:,:,Singlemap(reg45(i),5)));
end
for i = 1 : length(reg46)
    out.rmcor46(i) = corr2(ebc4.rm(:,:,Singlemap(reg46(i),4)),ebc6.rm(:,:,Singlemap(reg46(i),6)));
end
for i = 1 : length(reg47)
    out.rmcor47(i) = corr2(ebc4.rm(:,:,Singlemap(reg47(i),4)),ebc7.rm(:,:,Singlemap(reg47(i),7)));
end
for i = 1 : length(reg48)
    out.rmcor48(i) = corr2(ebc4.rm(:,:,Singlemap(reg48(i),4)),ebc8.rm(:,:,Singlemap(reg48(i),8)));
end
for i = 1 : length(reg49)
    out.rmcor49(i) = corr2(ebc4.rm(:,:,Singlemap(reg49(i),4)),ebc9.rm(:,:,Singlemap(reg49(i),9)));
end
for i = 1 : length(reg410)
    out.rmcor410(i) = corr2(ebc4.rm(:,:,Singlemap(reg410(i),4)),ebc10.rm(:,:,Singlemap(reg410(i),10)));
end


%Correlation sorting
%Keep tuning
for i = 1 : length(pass12)
    if ~isempty(pass12)
        out.keepRMcor12(i) = out.rmcor12(find(reg12 == pass12(i)));
        out.corindKeep12(i) = find(reg12 == pass12(i));
    end
end
for i = 1 : length(pass13)
    if ~isempty(pass13)
        out.keepRMcor13(i) = out.rmcor13(find(reg13 == pass13(i)));
        out.corindKeep13(i) = find(reg13 == pass13(i));
    end
end
for i = 1 : length(pass14)
    if ~isempty(pass14)
        out.keepRMcor14(i) = out.rmcor14(find(reg14 == pass14(i)));
        out.corindKeep14(i) = find(reg14 == pass14(i));
    end
end
for i = 1 : length(pass23)
    if ~isempty(pass23)
        out.keepRMcor23(i) = out.rmcor23(find(reg23 == pass23(i)));
        out.corindKeep23(i) = find(reg23 == pass23(i));
    end
end
for i = 1 : length(pass24)
    if ~isempty(pass24)
        out.keepRMcor24(i) = out.rmcor24(find(reg24 == pass24(i)));
        out.corindKeep24(i) = find(reg24 == pass24(i));
    end
end
for i = 1 : length(pass34)
    if ~isempty(pass34)
        out.keepRMcor34(i) = out.rmcor34(find(reg34 == pass34(i)));
        out.corindKeep34(i) = find(reg34 == pass34(i));
    end
end

%Lose tuning
c = 1;
for i = 1 : length(EBCind12)
    if isempty(find(pass12 == EBCind12(i)))
        out.lossRMcor12(c) = out.rmcor12(find(reg12 == EBCind12(i)));
        out.corind12(c) = find(reg12 == EBCind12(i));
        c = c+1;
    end
end
c = 1;
for i = 1 : length(EBCind13)
    if isempty(find(pass13 == EBCind13(i)))
        out.lossRMcor13(c) = out.rmcor13(find(reg13 == EBCind13(i)));
        out.corind13(c) = find(reg13 == EBCind13(i));
        c = c+1;
    end
end
c = 1;
for i = 1 : length(EBCind14)
    if isempty(find(pass14 == EBCind14(i)))
        out.lossRMcor14(c) = out.rmcor14(find(reg14 == EBCind14(i)));
        out.corind14(c) = find(reg14 == EBCind14(i));
        c = c+1;
    end
end
c = 1;
for i = 1 : length(EBCind23)
    if isempty(find(pass23 == EBCind23(i)))
        out.lossRMcor23(c) = out.rmcor23(find(reg23 == EBCind23(i)));
        out.corind23(c) = find(reg23 == EBCind23(i));
        c = c+1;
    end
end
c = 1;
for i = 1 : length(EBCind24)
    if isempty(find(pass24 == EBCind24(i)))
        out.lossRMcor24(c) = out.rmcor24(find(reg24 == EBCind24(i)));
        out.corind24(c) = find(reg24 == EBCind24(i));
        c = c+1;
    end
end
c = 1;
for i = 1 : length(EBCind43)
    if isempty(find(pass34 == EBCind43(i)))
        out.lossRMcor34(c) = out.rmcor34(find(reg34 == EBCind43(i)));
        out.corind43(c) = find(reg34 == EBCind43(i));
        c = c+1;
    end
end

%Gain tuning
c = 1;
for i = 1 : length(EBCind21)
    if isempty(find(pass12 == EBCind21(i)))
        out.gainRMcor12(c) = out.rmcor12(find(reg12 == EBCind21(i)));
        out.corind21(c) = find(reg12 == EBCind21(i));
        c = c+1;
    end
end
c = 1;
for i = 1 : length(EBCind31)
    if isempty(find(pass13 == EBCind31(i)))
        out.gainRMcor13(c) = out.rmcor13(find(reg13 == EBCind31(i)));
        out.corind31(c) = find(reg13 == EBCind31(i));
        c = c+1;
    end
end
c = 1;
for i = 1 : length(EBCind41)
    if isempty(find(pass14 == EBCind41(i)))
        out.gainRMcor14(c) = out.rmcor14(find(reg14 == EBCind41(i)));
        out.corind41(c) = find(reg14 == EBCind41(i));
        c = c+1;
    end
end
c = 1;
for i = 1 : length(EBCind32)
    if isempty(find(pass23 == EBCind32(i)))
        out.gainRMcor23(c) = out.rmcor23(find(reg23 == EBCind32(i)));
        out.corind32(c) = find(reg23 == EBCind32(i));
        c = c+1;
    end
end
c = 1;
for i = 1 : length(EBCind42)
    if isempty(find(pass24 == EBCind42(i)))
        out.gainRMcor24(c) = out.rmcor24(find(reg24 == EBCind42(i)));
        out.corind42(c) = find(reg24 == EBCind42(i));
        c = c+1;
    end
end
c = 1;
for i = 1 : length(EBCind34)
    if isempty(find(pass34 == EBCind34(i)))
        out.gainRMcor34(c) =out.rmcor34(find(reg34 == EBCind34(i)));
        out.corind34(c) = find(reg34 == EBCind34(i));
        c = c+1;
    end
end

out.EBCgain12 = length(EBCreg21) - length(pass12);
out.EBCloss12 = length(EBCreg12) - length(pass12);
out.EBCRegloss12 = length(fullpass1) - length(EBCreg12);
out.EBCRegloss21 = length(fullpass2) - length(EBCreg21);

out.EBCgain13 = length(EBCreg31) - length(pass13);
out.EBCloss13 = length(EBCreg13) - length(pass13);
out.EBCRegloss13 = length(fullpass1) - length(EBCreg13);
out.EBCRegloss31 = length(fullpass3) - length(EBCreg31);

out.EBCgain14 = length(EBCreg41) - length(pass14);
out.EBCloss14 = length(EBCreg14) - length(pass14);
out.EBCRegloss14 = length(fullpass1) - length(EBCreg14);
out.EBCRegloss41 = length(fullpass4) - length(EBCreg41);

out.EBCgain23 = length(EBCreg32) - length(pass23);
out.EBCloss23 = length(EBCreg23) - length(pass23);
out.EBCRegloss23 = length(fullpass2) - length(EBCreg23);
out.EBCRegloss32 = length(fullpass3) - length(EBCreg32);

out.EBCgain24 = length(EBCreg42) - length(pass24);
out.EBCloss24 = length(EBCreg24) - length(pass24);
out.EBCRegloss24 = length(fullpass2) - length(EBCreg24);
out.EBCRegloss42 = length(fullpass4) - length(EBCreg42);

out.EBCgain43 = length(EBCreg34) - length(pass34);
out.EBCloss43 = length(EBCreg43) - length(pass34);
out.EBCRegloss43 = length(fullpass4) - length(EBCreg43);
out.EBCRegloss34 = length(fullpass3) - length(EBCreg34);
%}