%% Define MRL and Firing rate cutoff for social target, identify cells which passed for left and right target and save renamed figures

function [] = passedSocialPreliminary_V2(p)

mrlThreshSep = 0.205;
firingThresh = 0.013;
mrlThreshSephalf = 0.41;
firingThreshHalf = 0.00655;
mrlThreshBoth = 0.25;
mrlThreshWall = 0.159;

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
        if ~isempty(find(strncmp(fnames,'Trial1Trial2ObjectData_D1A3_bothOb_V2',38),1))
            cd([folders{i} '\Trial1Trial2ObjectData_D1A3_bothOb_V2']);
            load('ms.mat')
            load('EBCstats.mat')
            subfolders = strsplit(folders{i},'\');
            for j = length(subfolders) : -1 : 1
                if ~isempty(find(contains(subfolders{j}, 'camkii-'),1)) && isempty(find(contains(subfolders{j}, 'trial'),1)) && isempty(find(contains(subfolders{j}, 'habituation'),1))
                    mouse = subfolders{j};
                    day = subfolders{j-1};
                    break
                end
            end
            
            %Seperate Targets
            passedSepL = Passed_Social_Criteria_V3(out.mrallOL,ms,mrlThreshSep,firingThresh);
            passedSepR = Passed_Social_Criteria_V3(out.mrallOR,ms,mrlThreshSep,firingThresh);
            
            passedSepL1 = Passed_Social_Criteria_V3(out.mrallO1L,ms,mrlThreshSephalf,firingThreshHalf);
            passedSepR1 = Passed_Social_Criteria_V3(out.mrallO1R,ms,mrlThreshSephalf,firingThreshHalf);
            
            passedSepL2 = Passed_Social_Criteria_V3(out.mrallO2L,ms,mrlThreshSephalf,firingThreshHalf);
            passedSepR2 = Passed_Social_Criteria_V3(out.mrallO2R,ms,mrlThreshSephalf,firingThreshHalf);
            
            %Both Targets
            passedBoth = Passed_Social_Criteria_V3(out.mrallO,ms,mrlThreshBoth,firingThresh);
            
            passedBoth1 = Passed_Social_Criteria_V3(out.mrallO1,ms,mrlThreshBoth,firingThreshHalf);
            passedBoth2 = Passed_Social_Criteria_V3(out.mrallO1,ms,mrlThreshBoth,firingThreshHalf);
            
            %Wall
            passedWall = Passed_Social_Criteria_V3(out.mrallW,ms,mrlThreshWall,firingThresh);
            
            passedWall1 = Passed_Social_Criteria_V3(out.mrallW1,ms,mrlThreshWall,firingThreshHalf);
            passedWall2 = Passed_Social_Criteria_V3(out.mrallW2,ms,mrlThreshWall,firingThreshHalf);
            
            %Trial1vsTrial2
            passedSepL_R = intersect(passedSepL,passedSepR);
            passedSepL1_L2 = intersect(passedSepL1,passedSepL2);
            passedSepL1_R1 = intersect(passedSepL1,passedSepR1);
            passedSepL1_R2 = intersect(passedSepL1,passedSepR2);
            passedSepL2_R1 = intersect(passedSepL2,passedSepR1);
            passedSepL2_R2 = intersect(passedSepL2,passedSepR2);
            passedSepR1_R2 = intersect(passedSepR1,passedSepR2);
            
            passedBoth1_2 = intersect(passedBoth1,passedBoth2);
            
            passedWall1_2 = intersect(passedWall1,passedWall2);
                     
            if ~isempty(passedWall1_2)
                for s = 1 : length(passedWall1_2)
                    copyfile(['./wall/wall' num2str(passedWall1_2(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\wall\Passed1&2\' day '_' mouse '_Cell' num2str(passedWall1_2(s)) '.jpg']);
                end
            end
            if ~isempty(passedWall1)
                for s = 1 : length(passedWall1)
                    copyfile(['./wall/wall' num2str(passedWall1(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\wall\Passed1\' day '_' mouse '_Cell' num2str(passedWall1(s)) '.jpg']);
                end
            end
            if ~isempty(passedWall2)
                for s = 1 : length(passedWall2)
                    copyfile(['./wall/wall' num2str(passedWall2(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\wall\Passed2\' day '_' mouse '_Cell' num2str(passedWall2(s)) '.jpg']);
                end
            end
            if ~isempty(passedWall)
                for s = 1 : length(passedWall)
                    copyfile(['./wall/wall' num2str(passedWall(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\wall\Passed\' day '_' mouse '_Cell' num2str(passedWall(s)) '.jpg']);
                end
            end
            
            if ~isempty(passedBoth1_2)
                for s = 1 : length(passedBoth1_2)
                    copyfile(['./ObjectsTogether/Objects' num2str(passedBoth1_2(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\BothTargets\Passed1&2\' day '_' mouse '_Cell' num2str(passedBoth1_2(s)) '.jpg']);
                end
            end
            if ~isempty(passedBoth1)
                for s = 1 : length(passedBoth1)
                    copyfile(['./ObjectsTogether/Objects' num2str(passedBoth1(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\BothTargets\Passed1\' day '_' mouse '_Cell' num2str(passedBoth1(s)) '.jpg']);
                end
            end
            if ~isempty(passedBoth2)
                for s = 1 : length(passedBoth2)
                    copyfile(['./ObjectsTogether/Objects' num2str(passedBoth2(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\BothTargets\Passed2\' day '_' mouse '_Cell' num2str(passedBoth2(s)) '.jpg']);
                end
            end
            if ~isempty(passedBoth)
                for s = 1 : length(passedBoth)
                    copyfile(['./ObjectsTogether/Objects' num2str(passedBoth(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\BothTargets\Passed\' day '_' mouse '_Cell' num2str(passedBoth(s)) '.jpg']);
                end
            end
            
            
            if ~isempty(passedSepL)
                for s = 1 : length(passedSepL)
                    copyfile(['./ObjectSeperate/Seperate' num2str(passedSepL(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\SeperateTargets\PassedLeft\' day '_' mouse '_Cell' num2str(passedSepL(s)) '.jpg']);
                end
            end
            %             end
            if ~isempty(passedSepR)
                for s = 1 : length(passedSepR)
                    copyfile(['./ObjectSeperate/Seperate' num2str(passedSepR(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\SeperateTargets\PassedRight\' day '_' mouse '_Cell' num2str(passedSepR(s)) '.jpg']);
                end
            end
            if ~isempty(passedSepL1)
                for s = 1 : length(passedSepL1)
                    copyfile(['./ObjectSeperate/Seperate' num2str(passedSepL1(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\SeperateTargets\PassedLeft1\' day '_' mouse '_Cell' num2str(passedSepL1(s)) '.jpg']);
                end
            end
            if ~isempty(passedSepL2)
                for s = 1 : length(passedSepL2)
                    copyfile(['./ObjectSeperate/Seperate' num2str(passedSepL2(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\SeperateTargets\PassedLeft2\' day '_' mouse '_Cell' num2str(passedSepL2(s)) '.jpg']);
                end
            end
            if ~isempty(passedSepR1)
                for s = 1 : length(passedSepR1)
                    copyfile(['./ObjectSeperate/Seperate' num2str(passedSepR1(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\SeperateTargets\PassedRight1\' day '_' mouse '_Cell' num2str(passedSepR1(s)) '.jpg']);
                end
            end
            if ~isempty(passedSepR2)
                for s = 1 : length(passedSepR2)
                    copyfile(['./ObjectSeperate/Seperate' num2str(passedSepR2(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\SeperateTargets\PassedRight2\' day '_' mouse '_Cell' num2str(passedSepR2(s)) '.jpg']);
                end
            end
            
            if ~isempty(passedSepL_R)
                for s = 1 : length(passedSepL_R)
                    copyfile(['./ObjectSeperate/Seperate' num2str(passedSepL_R(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\SeperateTargets\PassedLeftRight\' day '_' mouse '_Cell' num2str(passedSepL_R(s)) '.jpg']);
                end
            end
            
            if ~isempty(passedSepL1_R1)
                for s = 1 : length(passedSepL1_R1)
                    copyfile(['./ObjectSeperate/Seperate' num2str(passedSepL1_R1(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\SeperateTargets\PassedLeft1Right1\' day '_' mouse '_Cell' num2str(passedSepL1_R1(s)) '.jpg']);
                end
            end
            if ~isempty(passedSepL1_R2)
                for s = 1 : length(passedSepL1_R2)
                    copyfile(['./ObjectSeperate/Seperate' num2str(passedSepL1_R2(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\SeperateTargets\PassedLeft1Right2\' day '_' mouse '_Cell' num2str(passedSepL1_R2(s)) '.jpg']);
                end
            end
            if ~isempty(passedSepL2_R2)
                for s = 1 : length(passedSepL2_R2)
                    copyfile(['./ObjectSeperate/Seperate' num2str(passedSepL2_R2(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\SeperateTargets\PassedLeft2Right2\' day '_' mouse '_Cell' num2str(passedSepL2_R2(s)) '.jpg']);
                end
            end
            if ~isempty(passedSepL1_L2)
                for s = 1 : length(passedSepL1_L2)
                    copyfile(['./ObjectSeperate/Seperate' num2str(passedSepL1_L2(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\SeperateTargets\PassedLeft1Left2\' day '_' mouse '_Cell' num2str(passedSepL1_L2(s)) '.jpg']);
                end
            end
            if ~isempty(passedSepL2_R2)
                for s = 1 : length(passedSepL2_R2)
                    copyfile(['./ObjectSeperate/Seperate' num2str(passedSepL2_R2(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\SeperateTargets\PassedLeft2Right2\' day '_' mouse '_Cell' num2str(passedSepL2_R2(s)) '.jpg']);
                end
            end
            if ~isempty(passedSepL2_R1)
                for s = 1 : length(passedSepL2_R1)
                    copyfile(['./ObjectSeperate/Seperate' num2str(passedSepL2_R1(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\SeperateTargets\PassedLeft2Right1\' day '_' mouse '_Cell' num2str(passedSepL2_R1(s)) '.jpg']);
                end
            end
            if ~isempty(passedSepR1_R2)
                for s = 1 : length(passedSepR1_R2)
                    copyfile(['./ObjectSeperate/Seperate' num2str(passedSepR1_R2(s)) 'EBC.jpg'],['E:\SocialPrelimFigures\T1T2Comp_V2\SeperateTargets\PassedRight1Right2\' day '_' mouse '_Cell' num2str(passedSepR1_R2(s)) '.jpg']);
                end
            end
            save('PassedVals.mat', 'passedSepL','passedSepR','passedSepL1','passedSepR1','passedSepL2','passedSepR2','passedSepL_R','passedSepL1_R1','passedSepL1_R2','passedSepL2_R1','passedSepL2_R2','passedSepL1_L2','passedSepR1_R2','passedBoth','passedBoth1','passedBoth2','passedBoth1_2','passedWall','passedWall1','passedWall2','passedWall1_2')
        end
    end
end
end


