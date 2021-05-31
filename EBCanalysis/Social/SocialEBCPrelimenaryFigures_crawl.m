
function [] = SocialEBCPrelimenaryFigures_crawl(p)

paths = genpath(p);
if ispc
    folders = strsplit(paths,';')';
else
    folders = strsplit(paths,':')';
end
for i = 1 : length(folders)
    if ~isempty(folders{i})
        d = dir(folders{i});
        fnames = {d.name};
        if ~isempty(find(strcmp(fnames,'ms.mat'),1)) && ~isempty(find(strncmp(fnames,'EBCevenOddSplitParallelDeconvolved_Donut3SWEETMILK_D1A3_rightOb',1),1)) && ~isempty(find(strncmp(fnames,'EBCevenOddSplitParallelDeconvolved_Donut3SWEETMILK_D1A3.mat',1),1))
            cd(folders{i});                                  %Change current folder
            try
                load('ms.mat')
                load('EBCevenOddSplitParallelDeconvolved_Donut3SWEETMILK_D1A3.mat/EBCstats.mat')
                ebcl = out;
                load('EBCevenOddSplitParallelDeconvolved_Donut3SWEETMILK_D1A3_rightOb/EBCstats.mat')
                ebcr = out;
                passedSocialPreliminary(ebcl,ebcr,ms,0.41,0.00655);
            end
        end
    end
end
