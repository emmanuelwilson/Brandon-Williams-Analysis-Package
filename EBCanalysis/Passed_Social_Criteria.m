function [passed] = Passed_Social_Criteria(ebcoutput,calcium,mrlThresh,firingThresh)
% msses.deconvolvedSig(msses.deconvolvedSig > 0) = 1;
pMRL = find(ebcoutput.mrall>= mrlThresh);
if length(firingThresh) > 1
    pFiring = find(sum(calcium,1)/length(calcium(:,1))*5>=min(firingThresh) & sum(calcium,1)/length(calcium(:,1))*5 < max(firingThresh));
else 
    pFiring = find(sum(calcium,1)/length(calcium(:,1))*5>firingThresh);
end

passed = intersect(pMRL,pFiring);
end