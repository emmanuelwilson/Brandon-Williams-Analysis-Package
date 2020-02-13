%% Compares EBCs Ratemaps B context EBCs across A context EBCs
function [out] = EBCstabilityratemapcomp_AABA(Singlemap,fullpass1,fullpass2,fullpass3, fullpass4, ebc1, ebc2, ebc3, ebc4)
seshnum = 4;

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

EBCreg23 = Singlemap(ind2(find(Singlemap(ind2,3)>0)),3);
EBCreg24 = Singlemap(ind2(find(Singlemap(ind2,4)>0)),4);

EBCreg32 = Singlemap(ind3(find(Singlemap(ind3,2)>0)),2);
EBCreg34 = Singlemap(ind3(find(Singlemap(ind3,4)>0)),4);

EBCreg42 = Singlemap(ind4(find(Singlemap(ind4,2)>0)),2);
EBCreg43 = Singlemap(ind4(find(Singlemap(ind4,3)>0)),3);
EBCreg43 = Singlemap(ind4(find(Singlemap(ind4,3)>0)),3);


EBCind12 = ind1(find(Singlemap(ind1,2)>0));
EBCind21 = ind2(find(Singlemap(ind2,1)>0));
EBCind13 = ind1(find(Singlemap(ind1,3)>0));
EBCind31 = ind3(find(Singlemap(ind3,1)>0));
EBCind14 = ind1(find(Singlemap(ind1,4)>0));
EBCind41 = ind4(find(Singlemap(ind4,1)>0));

EBCind23 = ind2(find(Singlemap(ind2,3)>0));
EBCind24 = ind2(find(Singlemap(ind2,4)>0));

EBCind32 = ind3(find(Singlemap(ind3,2)>0));
EBCind34 = ind3(find(Singlemap(ind3,4)>0));

EBCind42 = ind4(find(Singlemap(ind4,2)>0));
EBCind43 = ind4(find(Singlemap(ind4,3)>0));

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
reg23 = find((bmap(:,2)+ bmap(:,3)) == 2);
reg24 = find((bmap(:,2)+ bmap(:,4)) == 2);
reg34 = find((bmap(:,4)+ bmap(:,3)) == 2);

CorrMap = cell(size(Singlemap,1),size(Singlemap,2),seshnum);
maxCorrMap = nan(size(Singlemap,1),size(Singlemap,2),seshnum);
medianCorrMap = nan(size(Singlemap,1),size(Singlemap,2),seshnum);
meanCorrMap = nan(size(Singlemap,1),size(Singlemap,2),seshnum);
minCorrMap = nan(size(Singlemap,1),size(Singlemap,2),seshnum);
maxrotCorrMap = nan(size(Singlemap,1),size(Singlemap,2),seshnum);
minrotCorrMap = nan(size(Singlemap,1),size(Singlemap,2),seshnum);

for mapind = 1 : length(Singlemap(:,1))
    if sum(bmap(mapind,:)) > 1
        for sessind1 = 1 : seshnum-1
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
                elseif sessind1 == 3
                    if Singlemap(mapind,4) > 0
                        CorrMap{mapind,sessind1,4} = rmRotcorr_singlecell(ebc3.rm(1:25,:,Singlemap(mapind,3)),ebc4.rm(:,:,Singlemap(mapind,4)));
                        [maxCorrMap(mapind,sessind1,4), maxrotCorrMap(mapind,sessind1,4)] = nanmax(CorrMap{mapind,sessind1,4});
                        medianCorrMap(mapind,sessind1,4) = nanmedian(CorrMap{mapind,sessind1,4});
                        meanCorrMap(mapind,sessind1,4) = nanmean(CorrMap{mapind,sessind1,4});
                        [minCorrMap(mapind,sessind1,4), minrotCorrMap(mapind,sessind1,4)] = nanmin(CorrMap{mapind,sessind1,4});
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

out.RegisteredMapIndex23 = reg23;
out.RegisteredMapIndex24 = reg24;

out.RegisteredMapIndex34 = reg34;


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