%% Compares the reliability of cells across days
% INPUT:
%   - folderpath: path to folder containing all pertinant information.
%   Namely:
%       *Singlemap: cell index across days 
%       *ShuffledCrit: reliability of cells and criteria structure
% OUTPUT: 
%   -out: ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Émmanuel Wilson

function out = ReliabilityAcrossDays_V2(folderpath)
oldpath = pwd;
cd(folderpath)
load('aligmentMaps','Singlemap')
Singlemap(Singlemap == 0) = NaN;
dCorrectRely = Singlemap;
dIncorrectRely = Singlemap;
dCorrectCorRely =Singlemap;
dIncorrectCorRely = Singlemap;
fCorrectRely = Singlemap;
fIncorrectRely = Singlemap;
fCorrectCorRely =Singlemap;
fIncorrectCorRely = Singlemap;
bCorrectRely = Singlemap;
bIncorrectRely = Singlemap;
bCorrectCorRely =Singlemap;
bIncorrectCorRely = Singlemap;

CritFiles = dir('Reliability/*ShuffledCrit*.mat');
for i = 1 : length(CritFiles)
    if isnumeric(str2num(CritFiles(i).name(13:end-4))) && (str2num(CritFiles(i).name(13:end-4)) == i)
        load(['Reliability/ShuffledCrit',num2str(i)],'ShuffledCrit')         
        for j = 1:length(ShuffledCrit.DelaySplithalfcorrect(:,1))
            dCorrectRely(find(Singlemap(:,i)== j),i) = ShuffledCrit.DelaySplithalfcorrect(j);
            dIncorrectRely(find(Singlemap(:,i)== j),i) = ShuffledCrit.DelaySplithalfincorrect(j);
            dCorrectCorRely(find(Singlemap(:,i)== j),i) = ShuffledCrit.DelaySplithalfccor(j);
            dIncorrectCorRely(find(Singlemap(:,i)== j),i) = ShuffledCrit.DelaySplithalficor(j);
            fCorrectRely(find(Singlemap(:,i)== j),i) = ShuffledCrit.FrontSplithalfcorrect(j);
            fIncorrectRely(find(Singlemap(:,i)== j),i) = ShuffledCrit.FrontSplithalfincorrect(j);
            fCorrectCorRely(find(Singlemap(:,i)== j),i) = ShuffledCrit.FrontSplithalfccor(j);
            fIncorrectCorRely(find(Singlemap(:,i)== j),i) = ShuffledCrit.FrontSplithalficor(j);
            bCorrectRely(find(Singlemap(:,i)== j),i) = ShuffledCrit.BackSplithalfcorrect(j);
            bIncorrectRely(find(Singlemap(:,i)== j),i) = ShuffledCrit.BackSplithalfincorrect(j);
            bCorrectCorRely(find(Singlemap(:,i)== j),i) = ShuffledCrit.BackSplithalfccor(j);
            bIncorrectCorRely(find(Singlemap(:,i)== j),i) = ShuffledCrit.BackSplithalficor(j);
        end               
    end
end

m = 5;
n = 4;
p = 1;
mkdir([folderpath,'/Reliability'],'Delay/CorrectReliabilityAcrossDays')
mkdir([folderpath,'/Reliability'],'Delay/IncorrectReliabilityAcrossDays')
mkdir([folderpath,'/Reliability'],'Delay/CorrectCorrectionReliabilityAcrossDays')
mkdir([folderpath,'/Reliability'],'Delay/IncorrectCorrectionReliabilityAcrossDays')

mkdir([folderpath,'/Reliability'],'FrontAnchored/CorrectReliabilityAcrossDays')
mkdir([folderpath,'/Reliability'],'FrontAnchored/IncorrectReliabilityAcrossDays')
mkdir([folderpath,'/Reliability'],'FrontAnchored/CorrectCorrectionReliabilityAcrossDays')
mkdir([folderpath,'/Reliability'],'FrontAnchored/IncorrectCorrectionReliabilityAcrossDays')

mkdir([folderpath,'/Reliability'],'BackAnchored/CorrectReliabilityAcrossDays')
mkdir([folderpath,'/Reliability'],'BackAnchored/IncorrectReliabilityAcrossDays')
mkdir([folderpath,'/Reliability'],'BackAnchored/CorrectCorrectionReliabilityAcrossDays')
mkdir([folderpath,'/Reliability'],'BackAnchored/IncorrectCorrectionReliabilityAcrossDays')

%% Delay
%Correct
for i = 1 : length(Singlemap(:,1))    
    if length(find(~isnan(dCorrectRely(i,:)))) > 1
        figure(1)
        subplot(m, n,p)
        p = p +1;
        plot(dCorrectRely(i,:))
        ylim([0 1])
        ylabel('Splithalf reliability')
        xlabel('Time (days)')
        title(['Correct Reliability Across Days, Cell: ',num2str(i)])
    end
    if p > 20 
        p = 1;
        saveas(gcf,['Reliability/Delay/CorrectReliabilityAcrossDays/CorrectRelyAcrossDays',num2str(i),'.jpg']);
        pause(0.01)
        clf
    end
end

%Incorrect
for i = 1 : length(Singlemap(:,1))    
    if length(find(~isnan(dIncorrectRely(i,:)))) > 1
        figure(1)
        subplot(m, n,p)
        p = p +1;
        plot(dIncorrectRely(i,:))
        ylim([0 1])
        ylabel('Splithalf reliability')
        xlabel('Time (days)')
        title(['Incorrect Reliability Across Days, Cell: ',num2str(i)])
    end
    if p > 20 
        p = 1;
        saveas(gcf,['Reliability/Delay/IncorrectReliabilityAcrossDays/IncorrectRelyAcrossDays',num2str(i),'.jpg']);
        pause(0.01)
        clf
    end
end

%Correct Correction
for i = 1 : length(Singlemap(:,1))    
    if length(find(~isnan(dCorrectCorRely(i,:)))) > 1
        figure(1)
        subplot(m, n,p)
        p = p +1;
        plot(dCorrectCorRely(i,:))
        ylim([0 1])
        ylabel('Splithalf reliability')
        xlabel('Time (days)')
        title(['Correct Correction Reliability Across Days, Cell: ',num2str(i)])
    end
    if p > 20 
        p = 1;
        saveas(gcf,['Reliability/Delay/CorrectCorrectionReliabilityAcrossDays/CorrectRelyAcrossDays',num2str(i),'.jpg']);
        pause(0.01)
        clf
    end
end

%Incorrect Correction
for i = 1 : length(Singlemap(:,1))    
    if length(find(~isnan(dIncorrectCorRely(i,:)))) > 1
        figure(1)
        subplot(m, n,p)
        p = p +1;
        plot(dIncorrectCorRely(i,:))
        ylim([0 1])
        ylabel('Splithalf reliability')
        xlabel('Time (days)')
        title(['Incorrect Correction Reliability Across Days, Cell: ',num2str(i)])
    end
    if p > 20 
        p = 1;
        saveas(gcf,['Reliability/Delay/IncorrectCorrectionReliabilityAcrossDays/CorrectRelyAcrossDays',num2str(i),'.jpg']);
        pause(0.01)
        clf
    end
end

%% Front Anchored
%Correct
for i = 1 : length(Singlemap(:,1))    
    if length(find(~isnan(fCorrectRely(i,:)))) > 1
        figure(1)
        subplot(m, n,p)
        p = p +1;
        plot(fCorrectRely(i,:))
        ylim([0 1])
        ylabel('Splithalf reliability')
        xlabel('Time (days)')
        title(['Correct Reliability Across Days, Cell: ',num2str(i)])
    end
    if p > 20 
        p = 1;
        saveas(gcf,['Reliability/FrontAnchored/CorrectReliabilityAcrossDays/CorrectRelyAcrossDays',num2str(i),'.jpg']);
        pause(0.01)
        clf
    end
end

%Incorrect
for i = 1 : length(Singlemap(:,1))    
    if length(find(~isnan(fIncorrectRely(i,:)))) > 1
        figure(1)
        subplot(m, n,p)
        p = p +1;
        plot(fIncorrectRely(i,:))
        ylim([0 1])
        ylabel('Splithalf reliability')
        xlabel('Time (days)')
        title(['Incorrect Reliability Across Days, Cell: ',num2str(i)])
    end
    if p > 20 
        p = 1;
        saveas(gcf,['Reliability/FrontAnchored/IncorrectReliabilityAcrossDays/IncorrectRelyAcrossDays',num2str(i),'.jpg']);
        pause(0.01)
        clf
    end
end

%Correct Correction
for i = 1 : length(Singlemap(:,1))    
    if length(find(~isnan(fCorrectCorRely(i,:)))) > 1
        figure(1)
        subplot(m, n,p)
        p = p +1;
        plot(fCorrectCorRely(i,:))
        ylim([0 1])
        ylabel('Splithalf reliability')
        xlabel('Time (days)')
        title(['Correct Correction Reliability Across Days, Cell: ',num2str(i)])
    end
    if p > 20 
        p = 1;
        saveas(gcf,['Reliability/FrontAnchored/CorrectCorrectionReliabilityAcrossDays/CorrectRelyAcrossDays',num2str(i),'.jpg']);
        pause(0.01)
        clf
    end
end

%Incorrect Correction
for i = 1 : length(Singlemap(:,1))    
    if length(find(~isnan(fIncorrectCorRely(i,:)))) > 1
        figure(1)
        subplot(m, n,p)
        p = p +1;
        plot(fIncorrectCorRely(i,:))
        ylim([0 1])
        ylabel('Splithalf reliability')
        xlabel('Time (days)')
        title(['Incorrect Correction Reliability Across Days, Cell: ',num2str(i)])
    end
    if p > 20 
        p = 1;
        saveas(gcf,['Reliability/FrontAnchored/IncorrectCorrectionReliabilityAcrossDays/CorrectRelyAcrossDays',num2str(i),'.jpg']);
        pause(0.01)
        clf
    end
end

%% Back Anchored
%Correct
for i = 1 : length(Singlemap(:,1))    
    if length(find(~isnan(bCorrectRely(i,:)))) > 1
        figure(1)
        subplot(m, n,p)
        p = p +1;
        plot(bCorrectRely(i,:))
        ylim([0 1])
        ylabel('Splithalf reliability')
        xlabel('Time (days)')
        title(['Correct Reliability Across Days, Cell: ',num2str(i)])
    end
    if p > 20 
        p = 1;
        saveas(gcf,['Reliability/BackAnchored/CorrectReliabilityAcrossDays/CorrectRelyAcrossDays',num2str(i),'.jpg']);
        pause(0.01)
        clf
    end
end

%Incorrect
for i = 1 : length(Singlemap(:,1))    
    if length(find(~isnan(bIncorrectRely(i,:)))) > 1
        figure(1)
        subplot(m, n,p)
        p = p +1;
        plot(bIncorrectRely(i,:))
        ylim([0 1])
        ylabel('Splithalf reliability')
        xlabel('Time (days)')
        title(['Incorrect Reliability Across Days, Cell: ',num2str(i)])
    end
    if p > 20 
        p = 1;
        saveas(gcf,['Reliability/BackAnchored/IncorrectReliabilityAcrossDays/IncorrectRelyAcrossDays',num2str(i),'.jpg']);
        pause(0.01)
        clf
    end
end

%Correct Correction
for i = 1 : length(Singlemap(:,1))    
    if length(find(~isnan(bCorrectCorRely(i,:)))) > 1
        figure(1)
        subplot(m, n,p)
        p = p +1;
        plot(bCorrectCorRely(i,:))
        ylim([0 1])
        ylabel('Splithalf reliability')
        xlabel('Time (days)')
        title(['Correct Correction Reliability Across Days, Cell: ',num2str(i)])
    end
    if p > 20 
        p = 1;
        saveas(gcf,['Reliability/BackAnchored/CorrectCorrectionReliabilityAcrossDays/CorrectRelyAcrossDays',num2str(i),'.jpg']);
        pause(0.01)
        clf
    end
end

%Incorrect Correction
for i = 1 : length(Singlemap(:,1))    
    if length(find(~isnan(bIncorrectCorRely(i,:)))) > 1
        figure(1)
        subplot(m, n,p)
        p = p +1;
        plot(bIncorrectCorRely(i,:))
        ylim([0 1])
        ylabel('Splithalf reliability')
        xlabel('Time (days)')
        title(['Incorrect Correction Reliability Across Days, Cell: ',num2str(i)])
    end
    if p > 20 
        p = 1;
        saveas(gcf,['Reliability/BackAnchored/IncorrectCorrectionReliabilityAcrossDays/CorrectRelyAcrossDays',num2str(i),'.jpg']);
        pause(0.01)
        clf
    end
end

%output
out.dCorrectRely = dCorrectRely;
out.dIncorrectRely = dIncorrectRely;
out.dCorrectCorRely = dCorrectCorRely;
out.dIncorrectCorRely = dIncorrectCorRely;
out.fCorrectRely = fCorrectRely;
out.fIncorrectRely = fIncorrectRely;
out.fCorrectCorRely = fCorrectCorRely;
out.fIncorrectCorRely = fIncorrectCorRely;
out.bCorrectRely = bCorrectRely;
out.bIncorrectRely = bIncorrectRely;
out.bCorrectCorRely = bCorrectCorRely;
out.bIncorrectCorRely = bIncorrectCorRely;
end