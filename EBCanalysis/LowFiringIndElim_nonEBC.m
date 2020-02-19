
function [nonEBCindnew] = LowFiringIndElim_nonEBC(ms1,ms2,ms3,ms4,firingT,fullpass1,fullpass2,fullpass3,fullpass4,nonEBCind,Singlemap)

fpass1 = find(sum(ms1.deconvolvedSig,1)/length(ms1.FiltTraces(:,1))*30>firingT);
fpass2 = find(sum(ms2.deconvolvedSig,1)/length(ms2.FiltTraces(:,1))*30>firingT);
fpass3 = find(sum(ms3.deconvolvedSig,1)/length(ms3.FiltTraces(:,1))*30>firingT);
fpass4 = find(sum(ms4.deconvolvedSig,1)/length(ms4.FiltTraces(:,1))*30>firingT);

[~,elim1,~] = intersect(fpass1,intersect(Singlemap(:,1), fullpass1));
[~,elim2,~] = intersect(fpass2,intersect(Singlemap(:,2), fullpass2));
[~,elim3,~] = intersect(fpass3,intersect(Singlemap(:,3), fullpass3));
[~,elim4,~] = intersect(fpass4,intersect(Singlemap(:,4), fullpass4));

fpass1(elim1) = [];
fpass2(elim2) = [];
fpass3(elim3) = [];
fpass4(elim4) = [];

[~,nonEBCind1,~] = intersect(Singlemap(nonEBCind,1),fpass1);
[~,nonEBCind2,~] = intersect(Singlemap(nonEBCind,2),fpass2);
[~,nonEBCind3,~] = intersect(Singlemap(nonEBCind,3),fpass3);
[~,nonEBCind4,~] = intersect(Singlemap(nonEBCind,4),fpass4);
nonEBCindnew = nonEBCind(unique([nonEBCind1;nonEBCind2;nonEBCind3;nonEBCind4]));

end
