%% Compares the mean Fluorescence of cells across days
% INPUT:
%   - folderpath: path to folder containing all pertinant information.
%   Namely:
%       *Singlemap: cell index across days
%       *Sep: Choice seperation of calcium for all trial types for specific
%           anchor point. 
% OUTPUT: 
%   -out: ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Émmanuel Wilson

function out = MeanFluoAcrossDays(folderpath,name,cut)
oldpath = pwd;
cd(folderpath)
load('aligmentMaps','Singlemap')
SepFiles = dir(['Reliability/*Sep',name,'*.mat']);
bmap = Singlemap;
bmap(bmap > 0) = 1;
elim = find(sum(bmap,2)==1);
ind = 1 : length(Singlemap(:,1));
ind(elim) = [];
max_est = 2000;
ChoiceNum = 5;
maxCor = 35;
count1 = cell(ChoiceNum,1);
count2 = cell(ChoiceNum,1);
count3 = cell(ChoiceNum,1);
count4 = cell(ChoiceNum,1);
CorrectComp = cell(ChoiceNum,1);
IncorrectComp = cell(ChoiceNum,1);
CorrectCorComp = cell(ChoiceNum,1);
IncorrectCorComp = cell(ChoiceNum,1);
behaveP = zeros(4,length(Singlemap(1,:)));
%All trials across days
CorrectAD = cell(ChoiceNum,1);
IncorrectAD = cell(ChoiceNum,1);
CorrectCorrectionAD = cell(ChoiceNum,1);
IncorrectCorrectionAD = cell(ChoiceNum,1);

for i = 1 : 5
    %Mean Fluorescence comparison across days
    CorrectComp{i} = zeros(length(ind),max_est,length(Singlemap(1,:)));
    IncorrectComp{i} = zeros(length(ind),max_est,length(Singlemap(1,:)));
    CorrectCorComp{i} = zeros(length(ind),max_est,length(Singlemap(1,:)));
    IncorrectCorComp{i} = zeros(length(ind),max_est,length(Singlemap(1,:)));
    
    %All trial comparison across days
    CorrectAD{i} = NaN(length(Singlemap(1,:))*maxCor,max_est,length(ind));
    IncorrectAD{i} = NaN(length(Singlemap(1,:))*maxCor,max_est,length(ind));
    CorrectCorrectionAD{i} = NaN(length(Singlemap(1,:))*maxCor,max_est,length(ind));
    IncorrectCorrectionAD{i} = NaN(length(Singlemap(1,:))*maxCor,max_est,length(ind));
end

for i = 1 : length(SepFiles)
    if (str2num(SepFiles(i).name((4+length(name)):end-4)) == i)
        sep =load(['Reliability/Sep',name,num2str(i)]);
        varname = fields(sep);
        sep = sep.(varname{1});
        subnames = fieldnames(sep);
        for j = 1 : length(ind)
            for c = 1 : ChoiceNum
                %Counter
                count1{c} = 1;
                count2{c} = 1;
                count3{c} = 1;
                count4{c} = 1;
            end
            cellnum = Singlemap(ind(j),i);
            if cellnum > 0
                for c = 1 : ChoiceNum
                    if any(strcmp(subnames, 'CorrectTrial')) && ~isempty(sep.CorrectTrial{c})
                        CorrectComp{c}(j,1:length(sep.CorrectTrial{c}(1,:,1)),i) = mean(sep.CorrectTrial{c}(:,:,cellnum),1);
                        CorrectAD{c}(count1{c} : count1{c} + length(sep.CorrectTrial{c}(:,1,1))-1,1:length(sep.CorrectTrial{c}(1,:,1)),j) = sep.CorrectTrial{c}(:,:,cellnum);
                        count1{c} = count1{c} + length(sep.CorrectTrial{c}(:,1,1));
                        if j ==1
                            behaveP(1,i) = behaveP(1,i) + length(sep.CorrectTrial{c}(:,1,1));
                        end
                    end
                    if any(strcmp(subnames, 'IncorrectTrial')) && ~isempty(sep.IncorrectTrial{c})
                        IncorrectComp{c}(j,1:length(sep.IncorrectTrial{c}(1,:,1)),i) = mean(sep.IncorrectTrial{c}(:,:,cellnum),1);
                        IncorrectAD{c}(count2{c} : count2{c} + length(sep.IncorrectTrial{c}(:,1,1))-1,1:length(sep.IncorrectTrial{c}(1,:,1)),j) = sep.IncorrectTrial{c}(:,:,cellnum);
                        count2{c} = count2{c} + length(sep.IncorrectTrial{c}(:,1,1));
                        if j == 1
                            behaveP(2,i) = behaveP(2,i) + length(sep.IncorrectTrial{c}(:,1,1));
                        end
                    end
                    if any(strcmp(subnames, 'CorrectCorrectionTrial')) && ~isempty(sep.CorrectCorrectionTrial{c})
                        CorrectCorComp{c}(j,1:length(sep.CorrectCorrectionTrial{c}(1,:,1)),i) = mean(sep.CorrectCorrectionTrial{c}(:,:,cellnum),1);
                        CorrectCorrectionAD{c}(count3{c} : count3{c} + length(sep.CorrectCorrectionTrial{c}(:,1,1))-1,1:length(sep.CorrectCorrectionTrial{c}(1,:,1)),j) = sep.CorrectCorrectionTrial{c}(:,:,cellnum);
                        count3{c} = count3{c} + length(sep.CorrectCorrectionTrial{c}(:,1,1));
                        if j == 1
                            behaveP(3,i) = behaveP(3,i) + length(sep.CorrectCorrectionTrial{c}(:,1,1));
                        end
                    end
                    if any(strcmp(subnames, 'IncorrectCorrectionTrial')) && ~isempty(sep.IncorrectCorrectionTrial{c})
                        IncorrectCorComp{c}(j,1:length(sep.IncorrectCorrectionTrial{c}(1,:,1)),i) = mean(sep.IncorrectCorrectionTrial{c}(:,:,cellnum),1);
                        IncorrectCorrectionAD{c}(count4{c} : count4{c} + length(sep.IncorrectCorrectionTrial{c}(:,1,1))-1,1:length(sep.IncorrectCorrectionTrial{c}(1,:,1)),j) = sep.IncorrectCorrectionTrial{c}(:,:,cellnum);
                        count4{c} = count4{c} + length(sep.IncorrectCorrectionTrial{c}(:,1,1));
                        if j == 1
                            behaveP(4,i) = behaveP(4,i) + length(sep.IncorrectCorrectionTrial{c}(:,1,1));
                        end
                    end
                end
            end
        end
    end
end

m1 = 0;
m2 = 0;

m = 5;
n = 4;

for i = 2 : ChoiceNum                
    if ~isempty(CorrectComp{i})
        if cut > 1
            mkdir([folderpath,'/Reliability/',name,'Anchored'],['MeanFluoAcrossDaysCorrect/',num2str(i)]);
            mkdir([folderpath,'/Reliability/',name,'Anchored'],['FluoAcrossDaysCorrect/',num2str(i)]);
        else
            mkdir([folderpath,'/Reliability/',name],['MeanFluoAcrossDaysCorrect/',num2str(i)]);
            mkdir([folderpath,'/Reliability/',name],['FluoAcrossDaysCorrect/',num2str(i)]);
        end
        for j = 1: length(CorrectComp{i}(1,1,:))
            m1 = length(find(CorrectComp{i}(1,:,j)));
            if m1 > m2
                m2 = m1;
            end
        end
        for j = length(CorrectAD{i}(:,1,1)) : -1 : 1
            if nansum(nansum(CorrectAD{i}(j,:,:))) == 0
                CorrectAD{i}(j,:,:) = [];
            end
        end
        CorrectComp{i} = CorrectComp{i}(:,1:m2,:);
        CorrectAD{i} = CorrectAD{i}(:,1:m2,:);
        if cut == 3 && length(CorrectComp{i}(1,:,1)) > 0
            for j = 1: length(CorrectComp{i}(1,1,:))
                for t = 1: length(CorrectComp{i}(:,1,1))
                    s = -(length(find(CorrectComp{i}(t,:,j))) - m2);
                    CorrectComp{i}(t,:,j) = circshift(CorrectComp{i}(t,:,j),s,2);
                end
            end
            if ~isempty(CorrectAD{i})
%                 for j = 1: length(CorrectAD{i}(1,1,:))
                    for t = 1: length(CorrectAD{i}(:,1,1))
                        s = -(length(find(~isnan(CorrectAD{i}(t,:,1)))) - m2);
                        CorrectAD{i}(t,:,:) = circshift(CorrectAD{i}(t,:,:),s,2);
                    end
%                 end
            end
        end
        c = 1;
        
        if any(find(CorrectComp{i}(:,:,:)))
            for j = 1 : length(ind)
                figure(1)
                subplot(m,n,c)
                c = c+1;
                im = permute(CorrectComp{i}(j,:,:),[3 2 1]);
                imagesc(im)
                colorbar
                title(['Cell ', num2str(ind(j))])
                if c > 20 || j == length(ind)
                    if cut > 1
                        saveas(gcf,[folderpath,'/Reliability/',name,'Anchored/MeanFluoAcrossDaysCorrect/',num2str(i),'/Cell', num2str(ind(j)),'CorrectTrials.jpg']);
                    else
                        saveas(gcf,[folderpath,'/Reliability/',name,'/MeanFluoAcrossDaysCorrect/',num2str(i),'/Cell', num2str(ind(j)),'CorrectTrials.jpg']);
                    end
                    pause(0.01)
                    clf
                    c = 1;
                end
            end
        end
        
        if any(find(CorrectAD{i}(:,:,:)))
            for j = 1 : length(ind)
                figure(1)
                subplot(m,n,c)
                c = c+1;               
                imagesc(CorrectAD{i}(:,:,j))
                colorbar
                title(['Cell ', num2str(ind(j))])
                if c > 20 || j == length(ind)
                    if cut > 1
                        saveas(gcf,[folderpath,'/Reliability/',name,'Anchored/FluoAcrossDaysCorrect/',num2str(i),'/Cell', num2str(ind(j)),'CorrectTrials.jpg']);
                    else
                        saveas(gcf,[folderpath,'/Reliability/',name,'/FluoAcrossDaysCorrect/',num2str(i),'/Cell', num2str(ind(j)),'CorrectTrials.jpg']);
                    end                    
                    pause(0.01)
                    clf
                    c = 1;
                end
            end
        end
            
    end
    %}
        m2 = 0;
    if ~isempty(IncorrectComp{i}) 
        if cut > 1
            mkdir([folderpath,'/Reliability/',name,'Anchored'],['MeanFluoAcrossDaysIncorrect/',num2str(i)]);
            mkdir([folderpath,'/Reliability/',name,'Anchored'],['FluoAcrossDaysIncorrect/',num2str(i)]);
        else
            mkdir([folderpath,'/Reliability/',name],['MeanFluoAcrossDaysIncorrect/',num2str(i)]);
            mkdir([folderpath,'/Reliability/',name],['FluoAcrossDaysIncorrect/',num2str(i)]);
        end
        for j = 1: length(IncorrectComp{i}(1,1,:))
            m1 = length(find(IncorrectComp{i}(1,:,j)));
            if m1 > m2
                m2 = m1;
            end
        end
        for j = length(IncorrectAD{i}(:,1,1)) : -1 : 1
            if nansum(nansum(IncorrectAD{i}(j,:,:))) == 0
                IncorrectAD{i}(j,:,:) = [];
            end
        end
        IncorrectComp{i} = IncorrectComp{i}(:,1:m2,:);
        IncorrectAD{i} = IncorrectAD{i}(:,1:m2,:);
        if cut == 3 && length(IncorrectComp{i}(1,:,1)) > 0
            for j = 1: length(IncorrectComp{i}(1,1,:))
                for t = 1 : length(IncorrectComp{i}(:,1,1))
                    s = -(length(find(IncorrectComp{i}(t,:,j))) - m2);
                    IncorrectComp{i}(t,:,j) = circshift(IncorrectComp{i}(t,:,j),s,2);
                end
            end
            if ~isempty(IncorrectAD{i})
%                 for j = 1: length(IncorrectAD{i}(1,1,:))
                    for t = 1 : length(IncorrectAD{i}(:,1,1))
                        s = -(length(find(~isnan(IncorrectAD{i}(t,:,1)))) - m2);
                        IncorrectAD{i}(t,:,:) = circshift(IncorrectAD{i}(t,:,:),s,2);
                    end
%                 end
            end
        end
        c = 1;
        if any(find(IncorrectComp{i}(:,:,:)))
            for j = 1 : length(ind)
                figure(1)
                subplot(m,n,c)
                c = c+1;
                im = permute(IncorrectComp{i}(j,:,:),[3 2 1]);
                imagesc(im)
                colorbar
                title(['Cell ', num2str(ind(j))])
                if c > 20 || j == length(ind)
                    if cut > 1
                        saveas(gcf,[folderpath,'/Reliability/',name,'Anchored/MeanFluoAcrossDaysIncorrect/',num2str(i),'/Cell', num2str(ind(j)),'IncorrectTrials.jpg']);
                    else
                        saveas(gcf,[folderpath,'/Reliability/',name,'/MeanFluoAcrossDaysIncorrect/',num2str(i),'/Cell', num2str(ind(j)),'IncorrectTrials.jpg']);
                    end                    
                    pause(0.01)
                    clf
                    c = 1;
                end
            end
        end
        if any(find(IncorrectAD{i}(:,:,:)))
            for j = 1 : length(ind)
                figure(1)
                subplot(m,n,c)
                c = c+1;                
                imagesc(IncorrectAD{i}(:,:,j))
                colorbar
                title(['Cell ', num2str(ind(j))])
                if c > 20 || j == length(ind)
                    if cut > 1
                        saveas(gcf,[folderpath,'/Reliability/',name,'Anchored/FluoAcrossDaysIncorrect/',num2str(i),'/Cell', num2str(ind(j)),'IncorrectTrials.jpg']);
                    else
                        saveas(gcf,[folderpath,'/Reliability/',name,'/FluoAcrossDaysIncorrect/',num2str(i),'/Cell', num2str(ind(j)),'IncorrectTrials.jpg']);
                    end                    
                    pause(0.01)
                    clf
                    c = 1;
                end
            end
        end
    end
    m2 = 0;
    if ~isempty(CorrectCorComp{i})
        if cut > 1
            mkdir([folderpath,'/Reliability/',name,'Anchored'],['MeanFluoAcrossDaysCorrectCorrection/',num2str(i)]);
            mkdir([folderpath,'/Reliability/',name,'Anchored'],['FluoAcrossDaysCorrectCorrection/',num2str(i)]);
        else
            mkdir([folderpath,'/Reliability/',name],['MeanFluoAcrossDaysCorrectCorrection/',num2str(i)]);
            mkdir([folderpath,'/Reliability/',name],['FluoAcrossDaysCorrectCorrection/',num2str(i)]);
        end
        for j = 1: length(CorrectCorComp{i}(1,1,:))
            m1 = length(find(CorrectCorComp{i}(1,:,j)));
            if m1 > m2
                m2 = m1;
            end
        end
        for j = length(CorrectCorrectionAD{i}(:,1,1)) : -1 : 1
            if nansum(nansum(CorrectCorrectionAD{i}(j,:,:))) == 0
                CorrectCorrectionAD{i}(j,:,:) = [];
            end
        end
        CorrectCorComp{i} = CorrectCorComp{i}(:,1:m2,:);
        CorrectCorrectionAD{i} = CorrectCorrectionAD{i}(:,1:m2,:);
        if cut == 3 && length(CorrectCorComp{i}(1,:,1)) > 0
            for j = 1: length(CorrectCorComp{i}(1,1,:))
                for t = 1: length(CorrectCorComp{i}(:,1,1))
                    s = -(length(find(CorrectCorComp{i}(t,:,j))) - m2);
                    CorrectCorComp{i}(t,:,j)= circshift(CorrectCorComp{i}(t,:,j),s,2);
                end
            end
            if ~isempty(CorrectCorrectionAD{i})
%                 for j = 1: length(CorrectCorrectionAD{i}(1,1,:))
                    for t = 1: length(CorrectCorrectionAD{i}(:,1,1))
                        s = -(length(find(~isnan(CorrectCorrectionAD{i}(t,:,1)))) - m2);
                        CorrectCorrectionAD{i}(t,:,:)= circshift(CorrectCorrectionAD{i}(t,:,:),s,2);
                    end
%                 end
            end
        end
        c = 1;
        if any(find(CorrectCorComp{i}(:,:,:)))
            for j = 1 : length(ind)
                figure(1)
                subplot(m,n,c)
                c = c+1;
                im = permute(CorrectCorComp{i}(j,:,:),[3 2 1]);
                imagesc(im)
                colorbar
                title(['Cell ', num2str(ind(j))])
                if c > 20 || j == length(ind)
                    if cut > 1
                        saveas(gcf,[folderpath,'/Reliability/',name,'Anchored/MeanFluoAcrossDaysCorrectCorrection/',num2str(i),'/Cell', num2str(ind(j)),'CorrectCorrectionTrials.jpg']);
                    else
                        saveas(gcf,[folderpath,'/Reliability/',name,'/MeanFluoAcrossDaysCorrectCorrection/',num2str(i),'/Cell', num2str(ind(j)),'CorrectCorrectionTrials.jpg']);
                    end                    
                    pause(0.01)
                    clf
                    c = 1;
                end
            end
        end
        if any(find(CorrectCorrectionAD{i}(:,:,:)))
            for j = 1 : length(ind)
                figure(1)
                subplot(m,n,c)
                c = c+1;                
                imagesc(CorrectCorrectionAD{i}(:,:,j))
                colorbar
                title(['Cell ', num2str(ind(j))])
                if c > 20 || j == length(ind)
                    if cut > 1
                        saveas(gcf,[folderpath,'/Reliability/',name,'Anchored/FluoAcrossDaysCorrectCorrection/',num2str(i),'/Cell', num2str(ind(j)),'CorrectCorrectionTrials.jpg']);
                    else
                        saveas(gcf,[folderpath,'/Reliability/',name,'/FluoAcrossDaysCorrectCorrection/',num2str(i),'/Cell', num2str(ind(j)),'CorrectCorrectionTrials.jpg']);
                    end                    
                    pause(0.01)
                    clf
                    c = 1;
                end
            end
        end
    end
    m2 = 0;
    if ~isempty(IncorrectCorComp{i})
        if cut > 1 
            mkdir([folderpath,'/Reliability/',name,'Anchored'],['MeanFluoAcrossDaysIncorrectCorrection/',num2str(i)]);    
            mkdir([folderpath,'/Reliability/',name,'Anchored'],['FluoAcrossDaysIncorrectCorrection/',num2str(i)]);    
        else
            mkdir([folderpath,'/Reliability/',name],['MeanFluoAcrossDaysIncorrectCorrection/',num2str(i)]);
            mkdir([folderpath,'/Reliability/',name],['FluoAcrossDaysIncorrectCorrection/',num2str(i)]);
        end
        for j = 1: length(IncorrectCorComp{i}(1,1,:))
            m1 = length(find(IncorrectCorComp{i}(1,:,j)));
            if m1 > m2
                m2 = m1;
            end
        end
        for j = length(IncorrectCorrectionAD{i}(:,1,1)) : -1 : 1
            if nansum(nansum(IncorrectCorrectionAD{i}(j,:,:))) == 0
                IncorrectCorrectionAD{i}(j,:,:) = [];
            end
        end
        IncorrectCorComp{i} = IncorrectCorComp{i}(:,1:m2,:);
        IncorrectCorrectionAD{i} = IncorrectCorrectionAD{i}(:,1:m2,:);
        if cut == 3 && length(IncorrectCorComp{i}(1,:,1)) > 0
            for j = 1: length(IncorrectCorComp{i}(1,1,:))
                for t = 1 : length(IncorrectCorComp{i}(:,1,1))                   
                    s = -(length(find(IncorrectCorComp{i}(t,:,j))) - m2);
                    IncorrectCorComp{i}(t,:,j) = circshift(IncorrectCorComp{i}(t,:,j),s,2);
                end
            end
            if ~isempty(IncorrectCorrectionAD{i})
%                 for j = 1: length(IncorrectCorrectionAD{i}(1,1,:))
                    for t = 1 : length(IncorrectCorrectionAD{i}(:,1,1))
                        s = -(length(all(find(IncorrectCorrectionAD{i}(t,:,1)))) - m2);
                        IncorrectCorrectionAD{i}(t,:,:) = circshift(IncorrectCorrectionAD{i}(t,:,:),s,2);
                    end
%                 end
            end
        end
        c = 1;
        if any(find(IncorrectCorComp{i}(:,:,:)))
            for j = 1 : length(ind)
                figure(1)
                subplot(m,n,c)
                c = c+1;
                im = permute(IncorrectCorComp{i}(j,:,:),[3 2 1]);
                imagesc(im)
                colorbar
                title(['Cell ', num2str(ind(j))])
                if c > 20 || j == length(ind)
                    if cut > 1
                        saveas(gcf,[folderpath,'/Reliability/',name,'Anchored/MeanFluoAcrossDaysIncorrectCorrection/',num2str(i),'/Cell', num2str(ind(j)),'IncorrectCorrectionTrials.jpg']);
                    else
                        saveas(gcf,[folderpath,'/Reliability/',name,'/MeanFluoAcrossDaysIncorrectCorrection/',num2str(i),'/Cell', num2str(ind(j)),'IncorrectCorrectionTrials.jpg']);
                    end                    
                    pause(0.01)
                    clf
                    c = 1;
                end
            end
        end
        if any(find(IncorrectCorrectionAD{i}(:,:,:)))
            for j = 1 : length(ind)
                figure(1)
                subplot(m,n,c)
                c = c+1;               
                imagesc(IncorrectCorrectionAD{i}(:,:,j))
                colorbar
                title(['Cell ', num2str(ind(j))])
                if c > 20 || j == length(ind)
                    if cut > 1
                        saveas(gcf,[folderpath,'/Reliability/',name,'Anchored/FluoAcrossDaysIncorrectCorrection/',num2str(i),'/Cell', num2str(ind(j)),'IncorrectCorrectionTrials.jpg']);
                    else
                        saveas(gcf,[folderpath,'/Reliability/',name,'/FluoAcrossDaysIncorrectCorrection/',num2str(i),'/Cell', num2str(ind(j)),'IncorrectCorrectionTrials.jpg']);
                    end                    
                    pause(0.01)
                    clf
                    c = 1;
                end
            end
        end
    end       
end

out.Correct = CorrectComp;
out.Incorrect = IncorrectComp;
out.CorrectCorrection = CorrectCorComp;
out.IncorrectCorrection = IncorrectCorComp;
out.BehavePerformance = behaveP;
end