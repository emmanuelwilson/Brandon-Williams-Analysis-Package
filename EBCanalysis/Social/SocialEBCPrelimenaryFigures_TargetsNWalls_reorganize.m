function [] = SocialEBCPrelimenaryFigures_TargetsNWalls_reorganize(p)

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
        if ~isempty(find(strncmp(fnames,'PassedL',1),1)) && ~isempty(find(strncmp(fnames,'PassedR.mat',1),1))
            cd(folders{i});
            load('PassedLeft/passedL.mat')
            load('PassedRight/passedR.mat')
            try
                load('PassedVals_clean.mat')
            catch
                load('PassedVals.mat')
            end
            passedboth = intersect(passedL,passedR);
            subfolders = strsplit(folders{i},'\');
            try
                for j = length(subfolders) : -1 : 1
                    if ~isempty(find(contains(subfolders{j}, 'camkii-'),1)) && isempty(find(contains(subfolders{j}, 'trial'),1)) && isempty(find(contains(subfolders{j}, 'habituation'),1))
                        newName = subfolders{j};
                        if isempty(find(strcmp(newName, mouse),1))                            
                            mouse = newName;
                            day = subfolders{j-1};
                        end
                        break
                    elseif ~isempty(find(contains(subfolders{j}, 'trial'),1))
                        if ~isempty(find(contains(subfolders{j}, 'trial1'),1))
                            tname = 'trial1';
                        else
                            tname = 'trial2';
                        end
                        if ~isempty(passedL)
                            for s = 1 : length(passedL)
                                copyfile(['../../Trial1Trial2ObjectData_D1A3_bothOb/ObjectSeperate/Seperate' num2str(passedL(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp\SeperateTargets\PassedLeft\' day mouse '_' tname '_Cell' num2str(passedL(s))]);
                            end
                            for s = 1 : length(passedR)
                                copyfile(['../../Trial1Trial2ObjectData_D1A3_bothOb/ObjectSeperate/Seperate' num2str(passedL(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp\SeperateTargets\PassedRight\' day mouse '_' tname '_Cell' num2str(passedR(s))]);
                            end
                            for s = 1 : length(passedboth)
                                copyfile(['../../Trial1Trial2ObjectData_D1A3_bothOb/ObjectSeperate/Seperate' num2str(passedL(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp\SeperateTargets\PassedBoth\' day mouse '_' tname '_Cell' num2str(passedboth(s))]);
                            end
                            for s = 1 : length(passed)
                                copyfile(['../../Trial1Trial2ObjectData_D1A3_bothOb/ObjectSeperate/Seperate' num2str(passedL(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp\SeperateTargets\PassedBoth\' day mouse '_' tname '_Cell' num2str(passedboth(s))]);
                                copyfile(['../../Trial1Trial2ObjectData_D1A3_bothOb/ObjectSeperate/Seperate' num2str(passedL(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp\BothTargets\PassedBoth\' day mouse '_' tname '_Cell' num2str(passedboth(s))]);
                            end
                                copyfile('PassedRight',['E:\SocialPrelimFigures\Trial' subfolders{j}(end) '/Right']);
                        elseif ~isempty(find(contains(subfolders{j}, 'habituation'),1))
                            break
                        end
                    end
                end
        
                for j = length(subfolders) : -1 : 1
                    if ~isempty(find(contains(subfolders{j}, 'trial'),1))
                        copyfile('PassedLeft',['E:\SocialPrelimFigures\Trial' subfolders{j}(end) '/Left']);
                        copyfile('PassedRight',['E:\SocialPrelimFigures\Trial' subfolders{j}(end) '/Right']);
                    elseif ~isempty(find(contains(subfolders{j}, 'habituation'),1))
                        copyfile('PassedLeft','E:\SocialPrelimFigures\Habituation/Left');
                        copyfile('PassedRight','E:\SocialPrelimFigures\Habituation/Right');
                    end
                end
            end
        end
    end
end
end

%{
        subfolders = strsplit(folders{i},'\');
        for j = length(subfolders) : -1 : 1
            if ~isempty(find(contains(subfolders{j}, 'camkii-'),1)) && isempty(find(contains(subfolders{j}, 'trial'),1)) && isempty(find(contains(subfolders{j}, 'habituation'),1))
                newName = subfolders{j};
                if isempty(find(strcmp(newName, mouse),1))
                    countM = countM +1;
                    mouse = newName;                                                            
                end                                
                break
            elseif ~isempty(find(contains(subfolders{j}, 'trial1'),1))
                load('PassedVals_clean.mat')
                
            elseif ~isempty(find(contains(subfolders{j}, 'habituation'),1))
                count = 1;
            end
        end
%}