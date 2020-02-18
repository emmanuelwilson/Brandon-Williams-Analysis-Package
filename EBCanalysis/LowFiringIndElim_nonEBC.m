
function [nonEBCindnew] = LowFiringIndElim_nonEBC(ms1,ms2,ms3,ms4,firingT,fullpass1,fullpass2,fullpass3,fullpass4,nonEBCind,Singlemap)

fpass21 = find(sum(ms1.deconvolvedSig,1)/length(ms1.FiltTraces(:,1))*30>firingT);
fpass22 = find(sum(ms2.deconvolvedSig,1)/length(ms2.FiltTraces(:,1))*30>firingT);
fpass23 = find(sum(ms3.deconvolvedSig,1)/length(ms3.FiltTraces(:,1))*30>firingT);
fpass24 = find(sum(ms4.deconvolvedSig,1)/length(ms4.FiltTraces(:,1))*30>firingT);

[~,elim1,~] = intersect(fpass21,intersect(Singlemap(:,1), fullpass1));
[~,elim2,~] = intersect(fpass22,intersect(Singlemap(:,2), fullpass2));
[~,elim3,~] = intersect(fpass23,intersect(Singlemap(:,3), fullpass3));
[~,elim4,~] = intersect(fpass24,intersect(Singlemap(:,4), fullpass4));

fpass21(elim1) = [];
fpass22(elim2) = [];
fpass23(elim3) = [];
fpass24(elim4) = [];

[~,nonEBCind21,~] = intersect(Singlemap(nonEBCind,1),fpass21);
[~,nonEBCind22,~] = intersect(Singlemap(nonEBCind,2),fpass22);
[~,nonEBCind23,~] = intersect(Singlemap(nonEBCind,3),fpass23);
[~,nonEBCind24,~] = intersect(Singlemap(nonEBCind,4),fpass24);
nonEBCindnew = unique([nonEBCind1,nonEBCind2,nonEBCind3,nonEBCind4]);

end
