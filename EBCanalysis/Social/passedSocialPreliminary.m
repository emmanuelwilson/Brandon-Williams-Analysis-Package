%% Define MRL and Firing rate cutoff for social target, identify cells which passed for left and right target and save renamed figures

function [] = passedSocialPreliminary(ebcLeft,ebcRight,ms,mrlThresh,firingThresh)

passedL = Passed_Social_Criteria_V2(ebcLeft,ms,mrlThresh,firingThresh);
passedR = Passed_Social_Criteria_V2(ebcRight,ms,mrlThresh,firingThresh);

nameR = 'EBCevenOddSplitParallelDeconvolved_Donut3SWEETMILK_D1A3_rightOb';
nameL = 'EBCevenOddSplitParallelDeconvolved_Donut3SWEETMILK_D1A3.mat';

folderR = dir(nameR);
folderL = dir(nameL);

path = strsplit(pwd, '\');
mousename = path{end-2};
date = path{end-3};

mkdir('PassedRight')
mkdir('PassedLeft')

for i = 1 : length(passedR)
    copyfile([nameR '/' num2str(passedR(i)) 'EBC.jpg'], ['PassedRight/' date mousename '_' num2str(passedR(i)) 'EBC.jpeg']);
end
for i = 1 : length(passedL)
    copyfile([nameL '/' num2str(passedL(i)) 'EBC.jpg'], ['PassedLeft/' date mousename '_' num2str(passedL(i)) 'EBC.jpeg']);
end

save('PassedLeft/passedL.mat', 'passedL')
save('PassedRight/passedR.mat', 'passedR')

end