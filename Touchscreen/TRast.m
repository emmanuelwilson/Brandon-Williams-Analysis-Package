%Creates raster plots of calcium activity for TouchScreen events marked with Scatterplots
%Specifically for: Front Anchored, Back Anchored and delay periods during
%Correct, Incorrect, Correct Correction and Incorrect Correction trials.
%
%INPUT: 
%- events: frameMap string file. Contains the synchronized frames, time
%   stamps of events as well as flags indicating event location.
%   events(:,1)= frameMap
%   events(:,2)= timestamp
%   events(:,3)= Correct/Incorrect/correction marker
%   events(:,5)= Trial start/finish markers
%   events(:,6)= Delay Markers
%-ms: Miniscope structure.
%   ms.FiltTraces = Filtered calcium traces of each cell
%-FA: logical 1 or 0 to either initiate or skip front anchored analysis
%-BA: logical 1 or 0 to either initiate or skip back anchored analysis
%-singleCells: logical 1 or 0 for single cell comparison accorss trials
%-splitd: 1 or 0 for split half comparison of delay periods.
%-savefigures: logical 1 or 0 to save the figures in the current folder. 

function TRast(events, ms, FA, BA,D, singleCells,SingleAllTrials,splitd, savefigures)
dirDelay = 'Delay';
dirFront = 'FrontAnchored';
dirBack = 'BackAnchored';


%seperate frame map from events file
for i = 1 : length(events)
    frameMap(i,1) = str2num(char(events(i,1)));
end

ctInd = [];                                                                 %Correct Trial Indices
itInd = [];                                                                 %Incorrect Trial Indices
ctCount = [];                                                               %Correct Trial count
itCount = [];                                                               %Incorrect Trial count
crcRow = [];                                                                %Correct correction trial row index
criRow = [];                                                                %Incorrect correction trial row index
crcInd = [];                                                                %Correct correction trial indices
criInd = [];                                                                %Incorrect correction trial indices
cri = 1;                                                                    %Correct correction trial index marker
crc = 1;                                                                    %Incorrect correction trial index marker
dc = [];
di = [];
dcc = [];
dic = [];
dccount = 1;
dicount = 1;
dcccount = 1;
diccount = 1;
Dcoef = zeros(length(ms.FiltTraces(1,:)),1);
Acoef = Dcoef;

%find all relevant events
Trial = find(contains(events(:,5),'1'));
Trial = [1 ; Trial];
correct = find(contains(events(:,3),'Correct'));
incorrect = find(contains(events(:,3),'Incorrect'));
correction = find(contains(events(:,3),'correction'));
delay = find(contains(events(:,6), '1'));
delay2 = find(contains(events(:,6), '2'));
timedelay = delay2-delay;
dFrames = mode(timedelay);
% delay(find(timedelay > (mode(timedelay)+1) || timedelay<(mode(timedelay)-1))) = [];
% delay2(find(timedelay > (mode(timedelay)+1) || timedelay<(mode(timedelay)-1))) = [];
% timedelay(find(timedelay > (mode(timedelay)+1) || timedelay<(mode(timedelay)-1))) = [];

%find nearest events to correct/incorrect choices
% [Cneighbr, ic,~] = unique(knnsearch(Trial,correct,'k',2),'rows');           %USE IC AND II TO CHANGE THE INDICES FOR THE SORTING BOUND LIMITS
Cneighbr = knnsearch(Trial,correct,'k',2);
Ineighbr = knnsearch(Trial,incorrect,'k',2);
% [Ineighbr, ii, ~] = unique(knnsearch(Trial,incorrect,'k',2),'rows');
Cneighbr = sort(Cneighbr,2);
Ineighbr = sort(Ineighbr,2);

flagInf = 0;

%Sort frames between correct/incorrect/correction trials
%Correct trial sorting
for i = 1 : length(Cneighbr)
    %Sorting lower bounds
    while min(Trial(Cneighbr(i,:)))> correct(i,1)
        temp = Cneighbr(find(Cneighbr == min(Cneighbr(i,:)),1)); %i == ic mod
        Cneighbr(find(Cneighbr == min(Cneighbr(i,1)),1)) = min(Cneighbr(i,:))-1;
        Cneighbr(find(Cneighbr(:,2) == max(Cneighbr(i,2)),1),2) = temp(1);
    end
    %Sorting upper bounds
    while ~any(isnan(Cneighbr(i,:)))&& flagInf <3
        if max(Trial(Cneighbr(i,:))) < correct(i,1)  %i == ic mod
            if(max(Cneighbr(i,:)))<= length(Trial)
                Cneighbr(i,find(max(Cneighbr(i,2)) == Cneighbr(i,:),1)) = NaN;
                Cneighbr(i,find(min(Cneighbr(i,1)) == Cneighbr(i,:),1)) = Cneighbr(i,find(min(Cneighbr(i,:)) == Cneighbr(i,1),1))+1;
            else
                Cneighbr(find(min(Cneighbr == Cneighbr(i,:)),1)) = Cneighbr(find(min(Cneighbr == Cneighbr(i,:)),1)+1);
            end
        end
        flagInf = flagInf +1;
    end
    ctInd(i,1) = Trial(nanmin(Cneighbr(i,:)));
%     if i>1 && any(ctInd(i,1) == ctInd(1:i-1,1))
%         ctInd(i,1) = 0;
%     end
    %Identifying Correct correction trials
    if any(ctInd(i,1) == correction(:,1))
        crcRow(i,1) = i;
%         crcInd(crc,1) = ctInd(i,1);
        crc = crc + 1;
    end
    %setting trial length
    if any(isnan(Cneighbr(i,:)))
        ctCount(i,1) = length(frameMap) - Trial(nanmax(Cneighbr(i,:)));
    else
        ctCount(i,1) = abs(diff(Trial(Cneighbr(i,:))));
    end
    flagInf = 0;
end

for i = 1: length(Ineighbr)
    %Sorting lower bounds
%     if i ==1
%         while min(Trial(Ineighbr(i,:)))> incorrect(i,1) %i == ii mod
%             temp = Ineighbr(find(Ineighbr == min(Ineighbr(i,:)),1));
%             Ineighbr(find(Ineighbr == min(Ineighbr(i,:)),1)) = min(Ineighbr(i,:))-1;
%             Ineighbr(find(Ineighbr == max(Ineighbr(i,:)),1)) = temp(1);
%         end
%         while ~any(isnan(Ineighbr(i,:)))&& flagInf <3
%             if max(Trial(Ineighbr(i,:))) < incorrect(i,1) %i == ii mod
%                 if(max(Ineighbr(i,:)))>= length(Trial)
%                     Ineighbr(i,find(max(Ineighbr(i,:)) == Ineighbr(i,:),1)) = NaN;
%                     Ineighbr(i,find(min(Ineighbr(i,:)) == Ineighbr(i,:),1)) = Ineighbr(i,find(min(Ineighbr(i,:)) == Ineighbr(i,:),1))+1;
%                 else
%                     Ineighbr(find(min(Ineighbr == Ineighbr(i,:)),1)) = Ineighbr(find(min(Ineighbr == Ineighbr(i,:)),1)+1);
%                 end
%             end
%             flagInf = flagInf +1;
%         end
%     else
        while min(Trial(Ineighbr(i,:)))> incorrect(i,1)% && ~(min(Ineighbr(i,:) == Ineighbr(i-1,:)))
            temp = Ineighbr(find(Ineighbr == min(Ineighbr(i,:)),1));            
            Ineighbr(find(Ineighbr == min(Ineighbr(i,:)),1)) = min(Ineighbr(i,:))-1;
            Ineighbr(find(Ineighbr(:,2) == max(Ineighbr(i,:)),1),2) = temp(1);
        end
        while ~any(isnan(Ineighbr(i,:)))&& flagInf <3
            if max(Trial(Ineighbr(i,:))) < incorrect(i,1)% && ~(min(Ineighbr(i,:) == Ineighbr(i-1,:)))
                if(max(Ineighbr(i,:)))>= length(Trial)
                    Ineighbr(i,find(max(Ineighbr(i,:)) == Ineighbr(i,:),1)) = NaN;
                    Ineighbr(i,find(min(Ineighbr(i,:)) == Ineighbr(i,:),1)) = Ineighbr(i,find(min(Ineighbr(i,:)) == Ineighbr(i,:),1))+1;
                else
                    Ineighbr(find(min(Ineighbr == Ineighbr(i,:)),1)) = Ineighbr(find(min(Ineighbr == Ineighbr(i,:)))+1);
                end
            end
            flagInf = flagInf +1;
        end
%     end
    
    itInd(i,1) = Trial(nanmin(Ineighbr(i,:)));
%     if i>1 && any(itInd(i,1) == itInd(1:i-1,1))
%         itInd(i,1) = 0;
%     end
    %Identifying Incorrect correction trials
    if any(itInd(i,1) == correction(:,1))
        criRow(i,1) = i;
%         criInd(cri,1) = itInd(i,:);
        cri = cri + 1;
    end
    %setting trial length
    if any(isnan(Ineighbr(i,:)))
        itCount(i,1) = length(frameMap) - Trial(nanmax(Ineighbr(i,:)));
    else
        itCount(i,1) = abs(diff(Trial(Ineighbr(i,:))));
    end
    flagInf = 0;
end
% 
[Cneighbr, ic,~] = unique(Cneighbr,'rows','stable');
[Ineighbr, ii, ~] = unique(Ineighbr,'rows','stable');
ctCount = ctCount(ic);%unique(ctCount,'rows','stable');
itCount = itCount(ii);%unique(itCount,'rows','stable');
ctInd = unique(ctInd,'stable');
itInd = unique(itInd,'stable');
% criInd = unique(criInd,'stable');
% crcInd = unique(crcInd,'stable');

%Slicing Calcium trace into their respective pieces,correct(Ctrace)/incorrect(Itrace)/correction(CrTrace)
% and taking out correction trials from Correct/Incorrect recordings
if ~isempty(crcRow)
%     crcRow = find(crcRow);
    [~,ia,~] = intersect(ic,crcRow);
    CrneighbrC = Cneighbr(ia,:);
    crCountC = ctCount(ia,:);
    crcInd = ctInd(ia);
    ctInd(ia,:)=  [];
    Cneighbr(ia,:)=  [];
    ctCount(ia,:)=  [];
end
if ~isempty(criRow)
%     criRow = find(criRow);
    [~,ia,~] = intersect(ii,criRow);
    CrneighbrI = Ineighbr(ia,:);
    crCountI = itCount(ia,:);    
    criInd = itInd(ia);
    itInd(ia,:)=  [];
    Ineighbr(ia,:)=  [];
    itCount(ia,:)=  [];
else
    CrneighbrI = zeros(1,2);
    crCountI = [];
end

%sort delay periods
sorted = 0;
for i = 1: length(delay)
    for j = 1 : length(Cneighbr(:,1))
        if delay(i) < Trial(max(Cneighbr(j,:))) && delay(i) > Trial(min(Cneighbr(j,:)))|| (any(find(isnan(Cneighbr(j,:)))) && delay(i) > Trial(min(Cneighbr(j,:))))
            dc(dccount) = delay(i);
            dccount = dccount +1;
            sorted = 1;
            break
        end
    end
    if ~sorted
        for j = 1 : length(Ineighbr(:,1))
            if delay(i) < Trial(max(Ineighbr(j,:))) && delay(i) > Trial(min(Ineighbr(j,:)))|| (any(find(isnan(Ineighbr(j,:)))) && delay(i) > Trial(min(Ineighbr(j,:))))
                di(dicount) = delay(i);
                dicount = dicount +1;
                sorted = 1;
                break
            end
        end
    end
    if ~sorted
        for j = 1 : length(CrneighbrC(:,1))
            if (delay(i) < Trial(max(CrneighbrC(j,:))) && delay(i) > Trial(min(CrneighbrC(j,:)))) || (any(find(isnan(CrneighbrC(j,:)))) && delay(i) > Trial(min(CrneighbrC(j,:))))
                dcc(dcccount) = delay(i);
                dcccount = dcccount +1;
                sorted = 1;
                break
            end
        end
    end
    if ~sorted && ~isempty(criRow)
        for j = 1 : length(CrneighbrI(:,1))
            if delay(i) < Trial(max(CrneighbrI(j,:))) && delay(i) > Trial(min(CrneighbrI(j,:))) || (any(find(isnan(CrneighbrI(j,:)))) && delay(i) > Trial(min(CrneighbrI(j,:))))
                dic(diccount) = delay(i);
                diccount = diccount +1;
                break
            end
        end
    end
    sorted = 0;
end

%Discecting the calcium traces in their respective trials
trace = transpose(normalize(ms.FiltTraces,'zscore'));
binarize = Binarize(ms);
btrace = transpose(binarize.binarizedTraces);
% trace = transpose(normalize(ms.FiltTraces,'range'));                        %total trace
Ctrace = nan(length(ms.FiltTraces(1,:)),max(ctCount),length(ctInd(:,1)));   %Correct Trace
Itrace = nan(length(ms.FiltTraces(1,:)),max(itCount),length(itInd(:,1)));   %Incorrect Trace
CRCtrace = nan(length(ms.FiltTraces(1,:)),max(crCountC),length(crcInd(:,1)));   %Correct correction trace
if ~isempty(criRow)
    CRItrace = nan(length(ms.FiltTraces(1,:)),max(crCountI),length(criInd(:,1)));   %Incorrect correction trace
end
dctrace = nan(length(ms.FiltTraces(1,:)),dFrames+1,length(dc));
ditrace = nan(length(ms.FiltTraces(1,:)),dFrames+1,length(di));
dcctrace = nan(length(ms.FiltTraces(1,:)),dFrames+1,length(dcc));
dicctrace = nan(length(ms.FiltTraces(1,:)),dFrames+1,length(dic));

for i = 1: length(ms.FiltTraces(1,:))
    for j =1: length(Ctrace(1,1,:))
        Ctrace(i,1:ctCount(j),j) = trace(i,ctInd(j):ctInd(j)+ctCount(j)-1);
    end    
    for j =1: length(Itrace(1,1,:))
        Itrace(i,1:itCount(j),j) = trace(i,itInd(j):itInd(j)+itCount(j)-1);
    end
    for j = 1: length(CRCtrace(1,1,:))
        CRCtrace(i,1:crCountC(j),j) = trace(i,crcInd(j):crcInd(j)+crCountC(j)-1);
    end
    if ~isempty(criRow)
        for j = 1: length(CRItrace(1,1,:))
            CRItrace(i,1:crCountI(j),j) = trace(i,criInd(j):criInd(j)+crCountI(j)-1);
        end
    else
        CRItrace = [];
    end
    for j =1: length(dc)
        dctrace(i,:,j) = trace(i,dc(j):dc(j)+dFrames); %Changed +60 to +dFrames
    end
    for j =1: length(di)
        ditrace(i,:,j) = trace(i,di(j):di(j)+dFrames);
    end
    for j =1: length(dcc)
        dcctrace(i,:,j) = trace(i,dcc(j):dcc(j)+dFrames);
    end
    for j =1: length(dic)
        dicctrace(i,:,j) = trace(i,dic(j):dic(j)+dFrames);
    end
    
    [Dcoef(i),Acoef(i)] = DelayActivityCoef(btrace(i,:),delay,dFrames);
end

% %Displaying individual trial delays for each cell
if singleCells && D   
%     if ~isempty(criRow)
    IndividualCells(dctrace,ditrace,dcctrace,dicctrace,ctCount,itCount,crCountC,crCountI,trace, 1,dirDelay, savefigures);   
%     else
%         IndividualCells(dctrace,ditrace,dcctrace,dicctrace,ctCount,itCount,crCountC,crCountI, 1,dirDelay, savefigures);   
%     end
end
if SingleAllTrials
%     salltrials = (dctrace,ditrace,dicctrace,ctCount,itCount,crcountC,crCountI,trace,,dirDelay,savefigures
end
if splitd
    delaySplitHalf(ms,dctrace,ditrace,dcctrace,dicctrace,dc,di,dcc,dic,delay,delay2);
end
% 
% %cut out excess/non-overlaping parts of the traces.
% if ~(min(ctCount) == length(Crast))
%     Crast = Crast(:,1:min(ctCount));
% end
% if ~(min(itCount) == length(Irast))
%     Irast = Irast(:,1:min(itCount));
% end
% if ~(min(crCountC) == length(CRCrast))
%     CRCrast = CRCrast(:,1:min(crCountC));
% end
% if ~(min(crCountI) == length(CRIrast))
%     CRIrast = CRIrast(:,1:min(crCountI));
% end
% 
% if singleCells
%     IndividualCells(Crast,Irast,CRCrast,CRIrast,dirFront);
% end
% 
% Crast = nanmean(Ctrace,3);
% Irast = nanmean(Itrace,3);
% CRCrast = nanmean(CRCtrace,3);
% CRIrast = nanmean(CRItrace,3);
if D
    dcrast = nanmean(dctrace,3);
    dirast = nanmean(ditrace,3);
    dccrast = nanmean(dcctrace,3);
    dicrast = nanmean(dicctrace,3);
    
    %Zero the minimum values Delay periods
    for i = 1 : length(dctrace(:,1,1))
        dcrast(i,:) = dcrast(i,:)-min(dcrast(i,:));
    end
    for i = 1 : length(ditrace(:,1,1))
        dirast(i,:) = dirast(i,:)-min(dirast(i,:));
    end
    for i = 1 : length(dcctrace(:,1,1))
        dccrast(i,:) = dccrast(i,:)-min(dccrast(i,:));
    end
    if ~isempty(criRow)
        for i = 1 : length(dicctrace(:,1,1))
            dicrast(i,:) = dicrast(i,:)-min(dicrast(i,:));
        end
    end
    
    dcrast = sortpeaks(dcrast);
    dirast = sortpeaks(dirast);
    dccrast = sortpeaks(dccrast);
    if ~isempty(criRow)
        dicrast = sortpeaks(dicrast);
    end
    
    %display the delay period
    d(1) = figure;
    h = imagesc(dcrast);
    % set(h,'LineStyle','none')
    axis('tight')
    colormap('hot')
    colorbar
    view(2)
    hold on
    xlabel('Time(frame)')
    ylabel('Neuron Number')
    title('Correct Trails Delay Period Calcium Activity')
    
    d(2) = figure;
    h = imagesc(dirast);
    % set(h,'LineStyle','none')
    axis('tight')
    colormap('hot')
    colorbar
    view(2)
    xlabel('Time(frame)')
    ylabel('Neuron Number')
    title('Incorrect Trails Delay Period Calcium Activity')
    
    d(3) = figure;
    h = imagesc(dccrast);
    % set(h,'LineStyle','none')
    axis('tight')
    colormap('hot')
    colorbar
    view(2)
    xlabel('Time(frame)')
    ylabel('Neuron Number')
    title('Correct Correction Trails Delay Period Calcium Activity')

    if ~isempty(criRow)
        d(4) = figure;
        h = imagesc(dicrast);
        % set(h,'LineStyle','none')
        axis('tight')
        colormap('hot')
        colorbar
        view(2)
        xlabel('Time(frame)')
        ylabel('Neuron Number')
        title('Incorrect Correction Trails Delay Period Calcium Activity')
    end
    
    if savefigures
        savefig(d,'Delay')
    end
end

% close(d)

if FA
    if singleCells
        IndividualCells(Ctrace,Itrace,CRCtrace,CRItrace,ctCount,itCount,crCountC,crCountI,trace,2,dirFront,savefigures);
    end
    Crast = nanmean(Ctrace,3);
    Irast = nanmean(Itrace,3);
    CRCrast = nanmean(CRCtrace,3);
    if ~isempty(criRow)
        CRIrast = nanmean(CRItrace,3);
    end
    
    %cut out excess/non-overlaping parts of the traces.
    %{    
    if ~(min(ctCount) == length(Crast))
        Crast = Crast(:,1:min(ctCount));
    end
    if ~(min(itCount) == length(Irast))
        Irast = Irast(:,1:min(itCount));
    end
    if ~(min(crCountC) == length(CRCrast))
        CRCrast = CRCrast(:,1:min(crCountC));
    end
    if ~isempty(criRow) && ~(min(crCountI) == length(CRIrast))
        CRIrast = CRIrast(:,1:min(crCountI));
    end
    %}
    
    %Zero the minimum values FrontAnchored
    for i = 1 : length(Ctrace(:,1))
        Crast(i,:) = Crast(i,:)-min(Crast(i,:));
    end
    for i = 1: length(Itrace(:,1))
        Irast(i,:) = Irast(i,:)-min(Irast(i,:));
    end
    if ~isempty(criRow)
        for i = 1: length(CRItrace(:,1))
            CRIrast(i,:) = CRIrast(i,:)-min(CRIrast(i,:));
        end
    end
    for i = 1: length(CRCtrace(:,1))
        CRCrast(i,:) = CRCrast(i,:)-min(CRCrast(i,:));
    end   
    
    %Sort cells by max peak time
    Crast = sortpeaks(Crast);
    Irast = sortpeaks(Irast);
    CRCrast = sortpeaks(CRCrast);
    if ~isempty(criRow)
        CRIrast = sortpeaks(CRIrast);
    end
    
    %display Front anchored averaging
    f(1) = figure;
    h = imagesc(Crast);
%     set(h,'LineStyle','none')
    axis('tight')
    colormap('hot')
    colorbar
    view(2)
    xlabel('Time(frame)')
    ylabel('Neuron Number')
    title('Correct Trails Calcium Activity, Front Anchored')
    
    f(2) = figure;
    h = imagesc(Irast);
%     set(h,'LineStyle','none')
    axis('tight')
    colormap('hot')
    colorbar
    view(2)
    title('Incorrect Trials Calcium Activity,Front Anchored')
    xlabel('Time(frame)')
    ylabel('Neuron Number')
    
    f(3) = figure;
    h = imagesc(CRCrast);
%     set(h,'LineStyle','none')
    axis('tight')
    colormap('hot')
    colorbar
    view(2)
    title('Correct Correction Trials Calcium Activity,Front Anchored')
    xlabel('Time(frame)')
    ylabel('Neuron Number')
    
    if ~isempty(criRow)
        f(4) = figure;
        h = imagesc(CRIrast);
%         set(h,'LineStyle','none')
        axis('tight')
        colormap('hot')
        colorbar
        view(2)
        title('Incorrect Correction Trials Calcium Activity,Front Anchored')
        xlabel('Time(frame)')
        ylabel('Neuron Number')
    end
    
    if savefigures
        savefig(f,'FrontAnchored')
    end
end
% close(f)

%reset traces and perform Back Anchored averaging
if BA
    Ctrace = nan(length(ms.FiltTraces(1,:)),max(ctCount),length(ctInd(:,1)));
    Itrace = nan(length(ms.FiltTraces(1,:)),max(itCount),length(itInd(:,1)));
    CRCtrace = nan(length(ms.FiltTraces(1,:)),max(crCountC),length(crcInd(:,1)));
    if ~isempty(criRow)
        CRItrace = nan(length(ms.FiltTraces(1,:)),max(crCountI),length(criInd(:,1)));
    end
    for i = 1: length(ms.FiltTraces(1,:))
        for j =1: length(Ctrace(1,1,:))
            Ctrace(i, (end- ctCount(j))+1:end,j) = trace(i,ctInd(j):ctInd(j)+ctCount(j)-1);
        end
        for j =1: length(Itrace(1,1,:))
            Itrace(i, (end- itCount(j))+1:end,j) = trace(i,itInd(j):itInd(j)+itCount(j)-1);
        end
        for j = 1: length(CRCtrace(1,1,:))
            CRCtrace(i,(end - crCountC(j))+1 :end,j) = trace(i,crcInd(j):crcInd(j)+crCountC(j)-1);
        end
        if ~isempty(criRow)
        for j = 1: length(CRItrace(1,1,:))
            CRItrace(i,(end-crCountI(j))+1 :end ,j) = trace(i,criInd(j):criInd(j)+crCountI(j)-1);
        end
        else 
            CRItrace = [];
        end
    end
    if singleCells
        IndividualCells(Ctrace,Itrace,CRCtrace,CRItrace,ctCount,itCount,crCountC,crCountI,trace,3, dirBack,savefigures);
    end
    
    Crast = nanmean(Ctrace,3);
    Irast = nanmean(Itrace,3);
    CRCrast = nanmean(CRCtrace,3);
    if ~isempty(criRow)
        CRIrast = nanmean(CRItrace,3);
    end
    
    %{
    if ~(min(ctCount) == length(Crast))
        Crast = Crast(:,end-min(ctCount):end);
    end       
    if ~(min(itCount) == length(Irast))
        Irast = Irast(:,end-min(itCount):end);
    end        
    if ~(min(crCountC) == length(CRCrast))
        CRCrast = CRCrast(:,end-min(crCountC):end);
    end        
    if ~isempty(criRow) && ~(min(crCountI) == length(CRIrast))
        CRIrast = CRIrast(:,end-min(crCountI):end);
    end
    %}  
    
    for i = 1 : length(Ctrace(:,1))
        Crast(i,:) = Crast(i,:)-min(Crast(i,:));
    end
    for i = 1: length(Itrace(:,1))
        Irast(i,:) = Irast(i,:)-min(Irast(i,:));
    end
    if ~isempty(criRow)
        for i = 1: length(CRItrace(:,1))
            CRIrast(i,:) = CRIrast(i,:)-min(CRIrast(i,:));
        end
    end
    for i = 1: length(CRCtrace(:,1))
        CRCrast(i,:) = CRCrast(i,:)-min(CRCrast(i,:));
    end
    
%     [~,maxind] = max(Crast,[],2);
%     [~, crastsort] = sort(maxind);
%     crastsort = flipud(crastsort);
%     Crast = Crast(crastsort,:);
%     [~,maxind] = max(Irast,[],2);
%     [~, irastsort] = sort(maxind);
%     irastsort = flipud(irastsort);
%     Irast = Irast(irastsort,:);
%     [~,maxind] = max(CRCrast,[],2);
%     [~, crcrastsort] = sort(maxind);
%     crcrastsort = flipud(crcrastsort);
%     CRCrast = CRCrast(crcrastsort,:);
%     [~,maxind] = max(CRIrast,[],2);
%     [~, crirastsort] = sort(maxind);
%     crirastsort = flipud(crirastsort);
%     CRIrast = CRIrast(crirastsort,:);
    
    Crast = sortpeaks(Crast);
    Irast = sortpeaks(Irast);
    CRCrast = sortpeaks(CRCrast);
    if ~isempty(criRow)
        CRIrast = sortpeaks(CRIrast);
    end
    
    %Display back anchored averaging raster plots
    b(1) = figure;
    h = imagesc(Crast);
%     set(h,'LineStyle','none')
    colormap('hot')
    colorbar
    axis('tight')
    view(2)
    xlabel('Time(frame)')
    ylabel('Neuron Number')
    title('Correct Trails Calcium Activity, Back Anchored')
    vline(length(Crast)-450,'g')
    
    b(2) = figure;
    h = imagesc(Irast);
%     set(h,'LineStyle','none')
    colormap('hot')
    colorbar
    axis('tight')
    view(2)
    title('Incorrect Trials Calcium Activity,Back Anchored')
    xlabel('Time(frame)')
    ylabel('Neuron Number')
    vline(length(Irast)-300,'b')
    
    b(3) = figure;
    h = imagesc(CRCrast);
%     set(h,'LineStyle','none')
    axis('tight')
    colormap('hot')
    colorbar
    view(2)
    title('Correct Correction Trials Calcium Activity,Back Anchored')
    xlabel('Time(frame)')
    ylabel('Neuron Number')
    vline(length(CRCrast)-450,'g')
    
    if ~isempty(criRow)
    b(4) = figure;
    h = imagesc(CRIrast);
%     set(h,'LineStyle','none')
    axis('tight')
    colormap('hot')
    colorbar
    view(2)
    title('Incorrect Correction Trials Calcium Activity,Back Anchored')
    xlabel('Time(frame)')
    ylabel('Neuron Number')
    vline(length(CRIrast)-300,'b')
    end
    if savefigures
        savefig(b,'BackAnchored')
    end
    % close(b)
end
end

%Will sort through the raster plot reoganizing cells by time of highest
%peak and return new raster plot
function rast = sortpeaks(rast)
[~,maxind] = max(rast,[],2);
[~, rastsort] = sort(maxind);
% rastsort = flipud(rastsort);
rast = rast(rastsort,:);
end

%This function will take your Correct Trial(ct), Incorrect Trial(it),
%Correct Correction trial(cct) and Incorrect Correction Trial matrices and
%compare each cell with itself across trials. It will save the results in
%respective subfolders within the folder indicated dirName.
function IndividualCells(ct,it,cct,ict,cl,il,ccl,icl,trace,cut,dirName,save)
figure(1)
n = 5;
m = 4;
minlim = min(trace,[],2);
maxlim = max(trace,[],2);
if ~isempty(ct) && length(ct(1,1,:))>1
    mkdir(dirName,'CorrectTrials');
    temp = permute(ct,[3 2 1]);
%     if cut == 1
%         %do nothing
%     elseif cut == 2
%         if ~(min(temp) == length(temp))
%             temp = temp(:,1:min(cl),:);
%         end
%     elseif cut == 3
%         if ~(min(temp) == length(temp))
%             temp = temp(:,end-min(cl):end,:);
%         end
%     end
    count = 1;
    %Correct trial comparison
    for i = 1 : length(ct(:,1))
        subplot(n,m,count)
        count = count +1;
        h = imagesc(temp(:,:,i));
%         set(h,'LineStyle','none')
        axis('tight')
        colormap('hot')
        colorbar
        caxis([minlim(i) maxlim(i)])
        view(2)
        xlabel('Time(frame)')
        ylabel('Trial Number')     
        title(['Cell ',num2str(i)])        
        hold on 
        for j = 1 : length(temp(:,1,1))
            if any(isnan(temp(j,:,i)))
                scatter(find(isnan(temp(j,:,i)),1,'first')-450,j,5,'o','g')
            else
                scatter(length((temp(j,:,i)))-450,j,5,'o','g')
            end
        end
        hold off
%         if cut == 3
%             vline(length(temp(1,:,1))-450,'g')
%         end
        if mod(i,n*m) == 0 || i == length(ct(:,1))            
            pause(0.01)
            count = 1;
            if save
                saveas(gcf,[dirName '/CorrectTrials/',num2str(i),'CorrectIndTrials.jpg']);
            end
        end
    end
end
if ~isempty(it) && length(it(1,1,:))>1
    mkdir(dirName,'IncorrectTrials');
    temp = permute(it,[3 2 1]);
%     if cut == 1
%         %do nothing
%     elseif cut == 2
%         if ~(min(temp) == length(temp))
%             temp = temp(:,1:min(il),:);
%         end
%     elseif cut == 3
%         if ~(min(temp) == length(temp))
%             temp = temp(:,end-min(il):end,:);
%         end
%     end
    count = 1;
    %Incorrect trial comparison
    for i = 1 : length(it(:,1))
        subplot(n,m,count)
        count = count +1;
        h = imagesc(temp(:,:,i));
%         set(h,'LineStyle','none')
        axis('tight')
        colormap('hot')
        colorbar
        caxis([minlim(i) maxlim(i)])
        view(2)
        xlabel('Time(frame)')
        ylabel('Trial Number')
        title(['Cell ',num2str(i)])
%         if cut == 3
%             vline(length(temp(1,:,1))-300,'b')
%         end
        hold on
        for j = 1 : length(temp(:,1,1))
            if any(isnan(temp(j,:,i)))
                scatter(find(isnan(temp(j,:,i)),1,'first')-300,j,5,'o','c')
            else
                scatter(length((temp(j,:,i)))-300,j,5,'o','c')
            end
        end
        hold off
        if mod(i,n*m) == 0 || i == length(it(:,1))
            pause(0.01)
            if save
                saveas(gcf,[dirName '/IncorrectTrials/',num2str(i),'IncorrectIndTrials.jpg']);
            end
            count = 1;
        end
    end
end
if ~isempty(cct) && length(cct(1,1,:))>1
    mkdir(dirName,'CorrectCorrection');
    temp = permute(cct,[3 2 1]);
%     if cut == 1
%         %do nothing
%     elseif cut == 2
%         if ~(min(temp) == length(temp))
%             temp = temp(:,1:min(ccl),:);
%         end
%     elseif cut == 3
%         if ~(min(temp) == length(temp))
%             temp = temp(:,end-min(ccl):end,:);
%         end
%     end
    count = 1;
    %Correct correction trial comparison
    for i = 1 : length(cct(:,1))
        subplot(n,m,count)
        count = count +1;
        h = imagesc(temp(:,:,i));
%         set(h,'LineStyle','none')
        axis('tight')
        colormap('hot')
        colorbar
        caxis([minlim(i) maxlim(i)])
        view(2)
        xlabel('Time(frame)')
        ylabel('Trial Number')
        title(['Cell ',num2str(i)])
        hold on
        for j = 1 : length(temp(:,1,1))
            if any(isnan(temp(j,:,i)))
                scatter(find(isnan(temp(j,:,i)),1,'first')-450,j,5,'o','g')
            else
                scatter(length((temp(j,:,i)))-450,j,5,'o','g')
            end
        end
        hold off
        if mod(i,n*m) == 0 || i == length(cct(:,1))
            pause(0.01)
            if save
                saveas(gcf,[dirName '/CorrectCorrection/',num2str(i),'CorrectCorrIndTrials.jpg']);
            end
            count = 1;
        end
    end
end
if ~isempty(ict) && length(ict(1,1,:))>1
    mkdir(dirName,'IncorrectCorrectTrials');
    temp = permute(ict,[3 2 1]);
%     if cut == 1
%         %do nothing
%     elseif cut == 2
%         if ~(min(temp) == length(temp))
%             temp = temp(:,1:min(icl),:);
%         end
%     elseif cut == 3
%         if ~(min(temp) == length(temp))
%             temp = temp(:,end-min(icl):end,:);
%         end
%     end
    count = 1;
    %Incorrect correction trial comparison
    for i = 1 : length(ict(:,1))        
        subplot(n,m,count)
        count = count +1;
        h = imagesc(temp(:,:,i));
%         set(h,'LineStyle','none')
        axis('tight')
        colormap('hot')
        colorbar
        caxis([minlim(i) maxlim(i)])
        view(2)
        xlabel('Time(frame)')
        ylabel('Trial Number')
        title(['Cell ',num2str(i)])
        
        hold on
        for j = 1 : length(temp(:,1,1))
            if any(isnan(temp(j,:,i)))
                scatter(find(isnan(temp(j,:,i)),1,'first')-450,j,5,'o','c')
            else
                scatter(length((temp(j,:,i)))-450,j,5,'o','c')
            end
        end
        hold off
        
        if mod(i,n*m) == 0 || i == length(ict(:,1))
            pause(0.01)
            if save 
                saveas(gcf,[dirName '/IncorrectCorrectTrials/',num2str(i),'IncorrectCorrIndTrials.jpg']);
            end
            count = 1;
        end
    end
end
end

function  activity = BinarizeRaster(traces)

% dt = median(diff(ms.time))/1000; % Conversion from ms to s
% Fs = 1/dt;
z_threshold = 2;

[bFilt,aFilt] = butter(2,  2/(30/2), 'low');

for trace_i = 1: length(traces(:,1))
    raw_trace = traces(trace_i,:);
    filt_trace = zscore(filtfilt(bFilt,aFilt,raw_trace));
    d1_trace = diff(filt_trace);
    d1_trace(end+1) = 0;
    d2_trace = diff(d1_trace);
    d2_trace(end+1) = 0;
    
    binary_trace = filt_trace*0;
    binary_trace(filt_trace>z_threshold & d1_trace>0) = 1;
    
    activity(trace_i,1) = sum(binary_trace,2);
    
    % spatialCoding.binarizedTraces(:, trace_i) = binary_trace;
    
end

end