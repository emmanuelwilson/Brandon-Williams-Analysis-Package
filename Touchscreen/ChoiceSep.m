%% Will sort trials based on mouse choice (ie.which square it clicked)
% INPUT: 
%   -events: Touchscreen sync variable
%       events is split as follows: 
%       events(:,1) = frameMap
%       events(:,2) = timestamp of event
%       events(:,3) = event name
%       events(:,4) = nose-poke position
%       events(:,5) = Trial Start/stop, indicated by a 1
%       events(:,6) = Delay period start/stop (1 = start, 2 = stop)
%   -cchoice: indices of correct choices (including correction trials)
%   -ichoice: indices of incorrect choices (including correction trials)
%   -cind: correct choice start of the trial index
%   -iind: incorrect choice start of the trial index
%   -ccind: correct correction choice start of the trial index
%   -icind: incorrect correction choice start of the trial index
%   -ct: correct trial traces
%   -it: incorrect trial traces
%   -cct: correct correction trial traces
%   -ict: incorrect correction trial traces
%   -cl: correct trial length 
%   -il: incorrect trial length 
%   -cl: correct correction trial length 
%   -icl: incorrect correction trial length 
%   -minmax: minimum and maxium calcium values for each cell
%   -cut: 1: Delay, 2: Front Anchored, 3: Back Anchored
%   -dirName: name of the folder in which the results will be saved in
%   -save: logical value, 1 save figures, 0 don't save figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Émmanuel Wilson

function [out,ccell,icell,cccell,iccell] = ChoiceSep(events, cchoice, ichoice, cind,iind,ccind,icind,ct,it,cct,ict,cl,il,ccl,icl,minmax,cut,dirName,save)
%% find which square was pressed during which trial type
cpress = zeros(length(cchoice),1);
ccpress = zeros(length(cchoice),1);
ipress = zeros(length(ichoice),1);
icpress = zeros(length(ichoice),1);

cc = 1;
ccc = 1;
ic = 1;
icc = 1;

for i = 1 : length(cchoice)    
    if length(cind)>= cc && (cchoice(i)>= cind(cc)) && (cchoice(i)<cind(cc)+cl(cc))
        cpress(i) = str2num(events(cchoice(i),4));
        cc = cc +1;
    elseif (cchoice(i)>= ccind(ccc)) && (cchoice(i)<ccind(ccc)+ccl(ccc))
        ccpress(i) = str2num(events(cchoice(i),4));
        ccc = ccc + 1;
    else
        warning(['missed choice', num2str(i)]);
    end
end
cpress = cpress(find(cpress));
ccpress = ccpress(find(ccpress));
for i = 1 : length(ichoice)
    if length(iind)>= ic && (ichoice(i)>= iind(ic)) && (ichoice(i)<iind(ic)+il(ic))
        ipress(i) = str2num(events(ichoice(i),4));
        ic = ic + 1;
    elseif ~isempty(icind) && (ichoice(i)>= icind(icc)) && (ichoice(i)<icind(icc)+icl(icc))
        icpress(i) = str2num(events(ichoice(i),4));
        icc = icc + 1; 
    else
        warning(['missed choice', num2str(i)]);
    end
end
ipress = ipress(find(ipress));
icpress = icpress(find(icpress));

%% sort which square was pressed during which trial
%sort choices and choice indices
[csort, ctind] = sort(cpress);
[ccsort, cctind] = sort(ccpress);
[isort, itind] = sort(ipress);
[icsort, ictind] = sort(icpress);
%Find how many different choices were made
cop = unique(cpress);
ccop = unique(ccpress);
iop = unique(ipress);
icop = unique(icpress);
%choice index cell matrix
ccell = cell(5,1);
cccell = ccell;
icell = ccell;
iccell = ccell;

ccount = 1;
cccount = 1;
icount = 1;
iccount = 1;
%fill cell array with corresponding choice trial indices
for i = 1 : 5
    if ccount <= length(cop) && cop(ccount) == i
        ind = find(csort == i);
        ccell{i} = ctind(ind);
        ccount = ccount + 1;
    end
    if cccount <= length(ccop) && ccop(cccount) == i
        ind = find(ccsort == i);
        cccell{i} = cctind(ind);
        cccount = cccount + 1;
    end
    if icount <= length(iop) && iop(icount) == i
        ind = find(isort == i);
        icell{i} = itind(ind);
        icount = icount + 1;
    end
    if ~isempty(icop) &&  iccount <= length(icop) && icop(iccount) == i
        ind = find(icsort == i);
        iccell{i} = ictind(ind);
        iccount = iccount + 1;
    end
end
%find total number of different choices
totpoke = unique(cat(1,cop,ccop,iop,icop));
%% Create matrices/graphs for sorted choice based trials
for posnum = 1 : length(totpoke)
    poke = totpoke(posnum);
    figure(1)
    n = 5;
    m = 4;
    minlim = minmax(:,1);
    maxlim = minmax(:,2);
    %Correct Trial
    if ~isempty(ct) &&  ~isempty(ccell{poke}) %%length(ct(1,1,:))>1 &&
        elim = 1:length(ct(1,1,:));
        elim(ccell{poke}) = [];
        mkdir(dirName,['ChoiceSeperation/CorrectTrials/', num2str(poke),'/AllCells']);
        temp = permute(ct,[3 2 1]);
        temp(elim,:,:) = [];
        if cut == 1
            %do nothing
        elseif cut == 2
            if ~(min(temp) == length(temp))
                temp = temp(:,1:min(cl),:);
            end
        elseif cut == 3
            if ~(min(temp) == length(temp))
                temp = temp(:,end+1-min(cl):end,:);
            end
        end
        out.CorrectTrial{poke} = temp;
        count = 1;
        if save%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Correct trial comparison
        for i = 1 : length(ct(:,1))
            subplot(n,m,count)
            count = count +1;
            h = imagesc(temp(:,:,i));
            %         set(h,'LineStyle','none')
            axis('tight')
            colormap('parula')
            colorbar
            caxis([minlim(i) maxlim(i)])
            view(2)
            xlabel('Time(frame)')
            ylabel('Trial Number')
            title(['Cell ',num2str(i)])
            if cut == 3
                vline(length(temp(1,:,1))-450,'g')
            end
            if mod(i,n*m) == 0 || i == length(ct(:,1))                
                count = 1;
                if save
                    saveas(gcf,[dirName '/ChoiceSeperation/CorrectTrials/', num2str(poke),'/AllCells/',num2str(i),'CorrectIndTrials.jpg']);
                    pause(0.01)
                end
                clf
            end
        end
        end%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    if isempty(ccell{poke})
        out.CorrectTrial{poke} = [];
    end
    %Incorrect Trial
    if ~isempty(it)  && ~isempty(icell{poke}) %%&& length(it(1,1,:))>1
        elim = 1:length(it(1,1,:));
        elim(icell{poke}) = [];
        mkdir(dirName,['ChoiceSeperation/IncorrectTrials/', num2str(poke),'/AllCells']);
        temp = permute(it,[3 2 1]);
        temp(elim,:,:) = [];
        if cut == 1
            %do nothing
        elseif cut == 2
            if ~(min(temp) == length(temp))
                temp = temp(:,1:min(il),:);
            end
        elseif cut == 3
            if ~(min(temp) == length(temp))
                temp = temp(:,end+1-min(il):end,:);
            end
        end
        out.IncorrectTrial{poke} = temp;
        count = 1;
        if save%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Incorrect trial comparison
        for i = 1 : length(it(:,1))
            subplot(n,m,count)
            count = count +1;
            h = imagesc(temp(:,:,i));
            %         set(h,'LineStyle','none')
            axis('tight')
            colormap('parula')
            colorbar
            caxis([minlim(i) maxlim(i)])
            view(2)
            xlabel('Time(frame)')
            ylabel('Trial Number')
            title(['Cell ',num2str(i)])
            if cut == 3
                vline(length(temp(1,:,1))-300,'b')
            end
            if mod(i,n*m) == 0 || i == length(it(:,1))                
                if save
                    saveas(gcf,[dirName '/ChoiceSeperation/IncorrectTrials/',num2str(poke),'/AllCells/',num2str(i),'IncorrectIndTrials.jpg']);
                    pause(0.01)
                end
                count = 1;
                clf
            end
        end
        end%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    if isempty(icell{poke})
        out.IncorrectTrial{poke} = [];
    end
    %Correct Correction
    if ~isempty(cct)  && ~isempty(cccell{poke})%%&& length(cct(1,1,:))>1
        elim = 1:length(cct(1,1,:));
        elim(cccell{poke}) = [];
        mkdir(dirName,['ChoiceSeperation/CorrectCorrectionTrials/', num2str(poke),'/AllCells']);
        temp = permute(cct,[3 2 1]);
        temp(elim,:,:) = [];
        if cut == 1
            %do nothing
        elseif cut == 2
            if ~(min(temp) == length(temp))
                temp = temp(:,1:min(ccl),:);
            end
        elseif cut == 3
            if ~(min(temp) == length(temp))
                temp = temp(:,end+1-min(ccl):end,:);
            end
        end
        out.CorrectCorrectionTrial{poke} = temp;
        count = 1;
        if save%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Correct correction trial comparison
        for i = 1 : length(cct(:,1))
            subplot(n,m,count)
            count = count +1;
            h = imagesc(temp(:,:,i));
            %         set(h,'LineStyle','none')
            axis('tight')
            colormap('parula')
            colorbar
            caxis([minlim(i) maxlim(i)])
            view(2)
            xlabel('Time(frame)')
            ylabel('Trial Number')
            title(['Cell ',num2str(i)])
            if cut == 3
                vline(length(temp(1,:,1))-450,'g')
            end
            if mod(i,n*m) == 0 || i == length(cct(:,1))                
                if save
                    saveas(gcf,[dirName '/ChoiceSeperation/CorrectCorrectionTrials/',num2str(poke),'/AllCells/',num2str(i),'CorrectCorrIndTrials.jpg']);
                    pause(0.01)
                end
                count = 1;
                clf
            end
        end
        end%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    if isempty(cccell{poke})
        out.CorrectCorrectionTrial{poke} = [];
    end
    %Incorrect Correction
    if ~isempty(ict)  && ~isempty(iccell{poke}) %&& length(ict(1,1,:))>1
        elim = 1:length(ict(1,1,:));
        elim(iccell{poke}) = [];
        mkdir(dirName,['ChoiceSeperation/IncorrectCorrectionTrials/', num2str(poke),'/AllCells']);
        temp = permute(ict,[3 2 1]);
        temp(elim,:,:) = [];
        if cut == 1
            %do nothing
        elseif cut == 2
            if ~(min(temp) == length(temp))
                temp = temp(:,1:min(icl),:);
            end
        elseif cut == 3
            if ~(min(temp) == length(temp))
                temp = temp(:,end+1-min(icl):end,:);
            end
        end
        out.IncorrectCorrectionTrial{poke} = temp;
        count = 1;
        if save%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Incorrect correction trial comparison
        for i = 1 : length(ict(:,1))
            subplot(n,m,count)
            count = count +1;
            h = imagesc(temp(:,:,i));
           %         set(h,'LineStyle','none')
            axis('tight')
            colormap('parula')
            colorbar
            caxis([minlim(i) maxlim(i)])
            view(2)
            xlabel('Time(frame)')
            ylabel('Trial Number')
            title(['Cell ',num2str(i)])
            if cut == 3
                vline(length(temp(1,:,1))-300,'b')
            end
            if mod(i,n*m) == 0 || i == length(ict(:,1))                
                if save
                    saveas(gcf,[dirName '/ChoiceSeperation/IncorrectCorrectionTrials/',num2str(poke),'/AllCells/',num2str(i),'IncorrectCorrIndTrials.jpg']);
                    pause(0.01)
                end
                count = 1;
                clf
            end
        end
        end%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    if isempty(iccell{poke})
        out.IncorrectCorrectionTrial{poke} = [];
    end
    
    %% Population level fluorescence 
    %Correct
    if length(out.CorrectTrial)>=poke && ~isempty(out.CorrectTrial{poke})
        raster = mean(out.CorrectTrial{poke},1);
        raster = permute(raster, [3 2 1]);
        [~,maxind] = max(raster,[],2);
        [~, rastsort] = sort(maxind);
        raster = raster(rastsort,:);
    else
        raster = [];
    end
    H(1) = figure;
    imagesc(raster)
    colormap parula
    colorbar
    if cut == 1
        title(['Delay: Correct Choice When selecting ', num2str(poke)])
    elseif cut == 2
        title(['Front Anchored: Correct Choice When selecting ', num2str(poke)])
    elseif cut == 3
        title(['Back Anchored: Correct Choice When selecting ', num2str(poke)])
    else
        title(['Correct Choice When selecting ', num2str(poke)])
    end
    xlabel('Time (frames)')
    ylabel('Cell')
    %Incorrect
    if  length(out.IncorrectTrial) >= poke && ~isempty(out.IncorrectTrial{poke})
        raster = mean(out.IncorrectTrial{poke},1);
        raster = permute(raster, [3 2 1]);
        [~,maxind] = max(raster,[],2);
        [~, rastsort] = sort(maxind);
        raster = raster(rastsort,:);        
    else
        raster = [];
    end
    H(2) = figure;
    imagesc(raster)
    colormap parula
    colorbar    
    if cut == 1
        title(['Delay: Incorrect Choice When selecting ', num2str(poke)])
    elseif cut == 2
        title(['Front Anchored: Incorrect Choice When selecting ', num2str(poke)])
    elseif cut == 3
        title(['Back Anchored: Incorrect Choice When selecting ', num2str(poke)])
    else
        title(['Incorrect Choice When selecting ', num2str(poke)])
    end
    xlabel('Time (frames)')
    ylabel('Cell')
    %Correct Correction 
    if length(out.CorrectCorrectionTrial) >= poke && ~isempty(out.CorrectCorrectionTrial{poke})
        raster = mean(out.CorrectCorrectionTrial{poke},1);
        raster = permute(raster, [3 2 1]);
        [~,maxind] = max(raster,[],2);
        [~, rastsort] = sort(maxind);
        raster = raster(rastsort,:);
    else
        raster = [];
    end
    H(3) = figure;
    imagesc(raster)
    colormap parula
    colorbar
    if cut == 1
        title(['Delay: Correct Correction Choice When selecting ', num2str(poke)])
    elseif cut == 2
        title(['Front Anchored: Correct Correction Choice When selecting ', num2str(poke)])
    elseif cut == 3
        title(['Back Anchored: Correct Correction Choice When selecting ', num2str(poke)])
    else
        title(['Correct Correction Choice When selecting ', num2str(poke)])
    end
    xlabel('Time (frames)')
    ylabel('Cell')
    %Incorrect Correction
    fnames = fieldnames(out);
    if any(strcmp(fnames, 'IncorrectCorrectionTrial')) && ~isempty(out.IncorrectCorrectionTrial{poke})
        raster = mean(out.IncorrectCorrectionTrial{poke},1);
        raster = permute(raster, [3 2 1]);
        [~,maxind] = max(raster,[],2);
        [~, rastsort] = sort(maxind);
        raster = raster(rastsort,:);
    else
        raster = [];
    end
    H(4) = figure;
    imagesc(raster)
    colormap parula
    colorbar   
    if cut == 1
        title(['Delay: Correct Incorrect Correction When selecting ', num2str(poke)])
    elseif cut == 2
        title(['Front Anchored: Incorrect Correction Choice When selecting ', num2str(poke)])
    elseif cut == 3
        title(['Back Anchored: Incorrect Correction Choice When selecting ', num2str(poke)])
    else
        title(['Incorrect Correction Choice When selecting ', num2str(poke)])
    end
    xlabel('Time (frames)')
    ylabel('Cell')
    if save 
        saveas(H,[dirName,'/Total_PopuplationChoice_',num2str(poke),'_AllCells'])
    end    
end




end