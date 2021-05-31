function [passed] = Passed_Social_Criteria_V2(ebcoutput,ms,mrlThresh,firingThresh)
% msses.deconvolvedSig(msses.deconvolvedSig > 0) = 1;
calcium = ms.deconvolvedSig;
if isfield(ms,'exclude')
    goodcell = find(ms.exclude.SFPs);
end
pMRL = find(ebcoutput.mrallO>= mrlThresh);
if length(firingThresh) > 1
    pFiring = find(sum(calcium,1)/length(calcium(:,1))*30>=min(firingThresh) & sum(calcium,1)/length(calcium(:,1))*30 < max(firingThresh));
else 
    pFiring = find(sum(calcium,1)/length(calcium(:,1))*30>firingThresh);
end

passed = intersect(pMRL,pFiring);
if isfield(ms,'exclude')
    passed = intersect(passed,goodcell);
end

end