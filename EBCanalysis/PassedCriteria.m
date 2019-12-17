% % % passed12 = find(out12.mrall>=0.1692);
% % % test = intersect(passed12,indOuttot12);
% % % for i  = length(test) : -1 : 1
% % % passed12(find(passed12==test(i))) = [];
% % % end

function [passed] = PassedCriteria(session,msses,mrlT,firingT)
% msses.deconvolvedSig(msses.deconvolvedSig > 0) = 1;
pMRL = find(session.mrall>= mrlT);
if length(firingT) > 1
    pFiring = find(sum(msses.deconvolvedSig,1)/length(msses.FiltTraces(:,1))*30>=min(firingT) & sum(msses.deconvolvedSig,1)/length(msses.FiltTraces(:,1))*30< max(firingT));
else 
    pFiring = find(sum(msses.deconvolvedSig,1)/length(msses.FiltTraces(:,1))*30>firingT);
end

passed = intersect(pMRL,pFiring);
end

% % test = intersect(passed21,indOuttot12);
% for i  = length(out21.mrall) : -1 : 1
% passed21(find(passed21==test(i))) = [];
% end
%{
passed101 = PassedCriteria(out101,ms101,0.1684,[0.0157, 0.03]);
passed102 = PassedCriteria(out102,ms102,0.1684,[0.0157, 0.03]);
passed103 = PassedCriteria(out103,ms103,0.1684,[0.0157, 0.03]);
passed104 = PassedCriteria(out104,ms104,0.1684,[0.0157, 0.03]);
passed106 = PassedCriteria(out106,ms106,0.1684,[0.0157, 0.03]);
passed108 = PassedCriteria(out108,ms108,0.1684,[0.0157, 0.03]);
passed1011 = PassedCriteria(out1011,ms1011,0.1684,[0.0157, 0.03]);
passed1012 = PassedCriteria(out1012,ms1012,0.1684,[0.0157, 0.03]);
passed1013 = PassedCriteria(out1013,ms1013,0.1684,[0.0157, 0.03]);
passed1014 = PassedCriteria(out1014,ms1014,0.1684,[0.0157, 0.03]);
passed24 = PassedCriteria(out24,ms24,0.1684,[0.0157, 0.03]);
passed25 = PassedCriteria(out25,ms25,0.1684,[0.0157, 0.03]);
passed26 = PassedCriteria(out26,ms26,0.1684,[0.0157, 0.03]);
passed28 = PassedCriteria(out28,ms28,0.1684,[0.0157, 0.03]);
passed2 = PassedCriteria(out29,ms29,0.1684,[0.0157, 0.03]);
passed29 = PassedCriteria(out29,ms29,0.1684,[0.0157, 0.03]);
passed210 = PassedCriteria(out210,ms210,0.1684,[0.0157, 0.03]);
passed211 = PassedCriteria(out211,ms211,0.1684,[0.0157, 0.03]);
passed212 = PassedCriteria(out212,ms212,0.1684,[0.0157, 0.03]);
passed214 = PassedCriteria(out214,ms214,0.1684,[0.0157, 0.03]);
passed21 = PassedCriteria(out21,ms21,0.1684,[0.0157, 0.03]);
passed22 = PassedCriteria(out22,ms22,0.1684,[0.0157, 0.03]);
passed23 = PassedCriteria(out23,ms23,0.1684,[0.0157, 0.03]);
passed61 = PassedCriteria(out61,ms61,0.1684,[0.0157, 0.03]);
passed62 = PassedCriteria(out62,ms62,0.1684,[0.0157, 0.03]);
passed63 = PassedCriteria(out63,ms63,0.1684,[0.0157, 0.03]);
passed64 = PassedCriteria(out64,ms64,0.1684,[0.0157, 0.03]);
passed65 = PassedCriteria(out65,ms65,0.1684,[0.0157, 0.03]);
passed69 = PassedCriteria(out69,ms69,0.1684,[0.0157, 0.03]);
passed610 = PassedCriteria(out610,ms610,0.1684,[0.0157, 0.03]);
passed611 = PassedCriteria(out611,ms611,0.1684,[0.0157, 0.03]);
passed612 = PassedCriteria(out612,ms612,0.1684,[0.0157, 0.03]);
passed613 = PassedCriteria(out613,ms613,0.1684,[0.0157, 0.03]);
passed614 = PassedCriteria(out614,ms614,0.1684,[0.0157, 0.03]);
passedA11 = PassedCriteria(outA11,msA11,0.1684,[0.0157, 0.03]);
passedA12 = PassedCriteria(outA12,msA12,0.1684,[0.0157, 0.03]);
passedA13 = PassedCriteria(outA13,msA13,0.1684,[0.0157, 0.03]);
passedA14 = PassedCriteria(outA14,msA14,0.1684,[0.0157, 0.03]);
passedA15 = PassedCriteria(outA15,msA15,0.1684,[0.0157, 0.03]);
passedA16 = PassedCriteria(outA16,msA16,0.1684,[0.0157, 0.03]);
passedA17 = PassedCriteria(outA17,msA17,0.1684,[0.0157, 0.03]);
passedA18 = PassedCriteria(outA18,msA18,0.1684,[0.0157, 0.03]);
passedA19 = PassedCriteria(outA19,msA19,0.1684,[0.0157, 0.03]);
passedA110 = PassedCriteria(outA110,msA110,0.1684,[0.0157, 0.03]);
passedA111 = PassedCriteria(outA111,msA111,0.1684,[0.0157, 0.03]);
passedA112 = PassedCriteria(outA112,msA112,0.1684,[0.0157, 0.03]);
passedA113 = PassedCriteria(outA113,msA113,0.1684,[0.0157, 0.03]);
passedA114 = PassedCriteria(outA114,msA114,0.1684,[0.0157, 0.03]);
passedA21 = PassedCriteria(outA21,msA21,0.1684,[0.0157, 0.03]);
passedA22 = PassedCriteria(outA22,msA22,0.1684,[0.0157, 0.03]);
passedA24 = PassedCriteria(outA24,msA24,0.1684,[0.0157, 0.03]);
passedA25 = PassedCriteria(outA25,msA25,0.1684,[0.0157, 0.03]);
passedA27 = PassedCriteria(outA27,msA27,0.1684,[0.0157, 0.03]);
passedA28 = PassedCriteria(outA28,msA28,0.1684,[0.0157, 0.03]);
passedA29 = PassedCriteria(outA29,msA29,0.1684,[0.0157, 0.03]);
passedA210 = PassedCriteria(outA210,msA210,0.1684,[0.0157, 0.03]);
passedA211 = PassedCriteria(outA211,msA211,0.1684,[0.0157, 0.03]);
passedA212 = PassedCriteria(outA212,msA212,0.1684,[0.0157, 0.03]);
passedA213 = PassedCriteria(outA213,msA213,0.1684,[0.0157, 0.03]);
passedA214 = PassedCriteria(outA214,msA214,0.1684,[0.0157, 0.03]);
passedA21 = PassedCriteria(outA21,msA21,0.1684,[0.0157, 0.03]);
passedA22 = PassedCriteria(outA22,msA22,0.1684,[0.0157, 0.03]);
passedA24 = PassedCriteria(outA24,msA24,0.1684,[0.0157, 0.03]);
passedA25 = PassedCriteria(outA25,msA25,0.1684,[0.0157, 0.03]);
passedA27 = PassedCriteria(outA27,msA27,0.1684,[0.0157, 0.03]);
passedA28 = PassedCriteria(outA28,msA28,0.1684,[0.0157, 0.03]);
passedA29 = PassedCriteria(outA29,msA29,0.1684,[0.0157, 0.03]);
passedA210 = PassedCriteria(outA210,msA210,0.1684,[0.0157, 0.03]);
passedA211 = PassedCriteria(outA211,msA211,0.1684,[0.0157, 0.03]);
passedA212 = PassedCriteria(outA212,msA212,0.1684,[0.0157, 0.03]);
passedA213 = PassedCriteria(outA213,msA213,0.1684,[0.0157, 0.03]);
passedA214 = PassedCriteria(outA214,msA214,0.1684,[0.0157, 0.03]);
passedW1 = PassedCriteria(outW1,msW1,0.1684,[0.0157, 0.03]);
passedW2 = PassedCriteria(outW2,msW2,0.1684,[0.0157, 0.03]);
passedW3 = PassedCriteria(outW3,msW3,0.1684,[0.0157, 0.03]);
passedW4 = PassedCriteria(outW4,msW4,0.1684,[0.0157, 0.03]);
passedW5 = PassedCriteria(outW5,msW5,0.1684,[0.0157, 0.03]);
passedW6 = PassedCriteria(outW6,msW6,0.1684,[0.0157, 0.03]);
passedW7 = PassedCriteria(outW7,msW7,0.1684,[0.0157, 0.03]);
passedW8 = PassedCriteria(outW8,msW8,0.1684,[0.0157, 0.03]);
passedW9 = PassedCriteria(outW9,msW9,0.1684,[0.0157, 0.03]);
passedW14 = PassedCriteria(outW14,msW14,0.1684,[0.0157, 0.03]);
passed805 = PassedCriteria(out805,ms805,0.1684,[0.0157, 0.03]);
passed806 = PassedCriteria(out805,ms806,0.1684,[0.0157, 0.03]);
passed806 = PassedCriteria(out806,ms806,0.1684,[0.0157, 0.03]);
passed807 = PassedCriteria(out807,ms807,0.1684,[0.0157, 0.03]);
%}