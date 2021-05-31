function [out] = passedSocialPreliminary_V2_Stats(p)

paths = genpath(p);
if ispc
    folders = strsplit(paths,';')';
else
    folders = strsplit(paths,':')';
end

out.Sep_MouseLeft2RightFollow = [];
out.Sep_MouseRight2LeftFollow = [];
out.Sep_CupLeftFollow = [];
out.Sep_CupRightFollow = [];
count = 1;

for i = 1 : length(folders)
    if ~isempty(folders{i})
        d = dir(folders{i});
        fnames = {d.name};
        if ~isempty(find(strncmp(fnames,'Trial1Trial2ObjectData_D1A3_bothOb',1),1))
            cd([folders{i} '\Trial1Trial2ObjectData_D1A3_bothOb']);
            load('PassedVals')
            
            
            out.Left1passnum(count) = length(passedSepL1);
            out.Left2passnum(count) = length(passedSepL2);
            
            out.Right1passnum(count) = length(passedSepR1);
            out.Right2passnum(count) = length(passedSepR2);
            
            out.Sep_MouseLeft2RightFollow(count) = length(passedSepL1_R2);
            out.Sep_MouseRight2LeftFollow(count) = length(passedSepL2_R1);
            out.Sep_CupLeftFollow(count) = length(passedSepL1_L2);
            out.Sep_CupRightFollow(count) = length(passedSepR1_R2);
            
            count = count + 1;
        end
    end
end