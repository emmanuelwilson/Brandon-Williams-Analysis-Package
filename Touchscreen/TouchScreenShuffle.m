%Shuffle Touchscreen data 100 times and compares to split half correlation
%
%INPUT:
%   -ms: miniscope structure, pertinant variables include:
%       *FiltTraces: n X m matrix, where n is the frame number and m the
%       cell number.
%   -events: string containing touchscreen event information
%OUTPUT:
%   - Mean shuffled splithalf correlation, mean random trial interval
%   split half, 99th percentile, # of cells passed, cell index of passed
%   cells, calcium of cells which passed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%By: Émmanuel Wilson 

function [out] = TouchScreenShuffle(ms,events,calInfo)
%seperate frame map from events file
for i = 1 : length(events)
    frameMap(i,1) = str2num(char(events(i,1)));
end
tic
%Initilize variables
out.ScorrDelaycorrect = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrDelayincorrect = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrDelayccor = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrDelayicor = zeros(length(ms.FiltTraces(1,:)),100);

out.DelaySplithalfcorrect = zeros(length(ms.FiltTraces(1,:)),1);
out.DelaySplithalfincorrect = zeros(length(ms.FiltTraces(1,:)),1);
out.DelaySplithalfccor = zeros(length(ms.FiltTraces(1,:)),1);
out.DelaySplithalficor = zeros(length(ms.FiltTraces(1,:)),1);

out.ScorrFrontcorrect = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrFrontincorrect = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrFrontccor = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrFronticor = zeros(length(ms.FiltTraces(1,:)),100);

out.FrontSplithalfcorrect = zeros(length(ms.FiltTraces(1,:)),1);
out.FrontSplithalfincorrect = zeros(length(ms.FiltTraces(1,:)),1);
out.FrontSplithalfccor = zeros(length(ms.FiltTraces(1,:)),1);
out.FrontSplithalficor = zeros(length(ms.FiltTraces(1,:)),1);

out.ScorrBackcorrect = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrBackincorrect = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrBackccor = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrBackicor = zeros(length(ms.FiltTraces(1,:)),100);

out.BackSplithalfcorrect = zeros(length(ms.FiltTraces(1,:)),1);
out.BackSplithalfincorrect = zeros(length(ms.FiltTraces(1,:)),1);
out.BackSplithalfccor = zeros(length(ms.FiltTraces(1,:)),1);
out.BackSplithalficor = zeros(length(ms.FiltTraces(1,:)),1);

ctInd = calInfo.correctStartInd;                                                                 %Correct Trial Indices
itInd = calInfo.incorrectStartInd;                                                                 %Incorrect Trial Indices
crcInd = calInfo.correctcorStartInd;                                                                %Correct correction trial indices
criInd = calInfo.incorrectcorStartInd;                                                                %Incorrect correction trial indices
ctCount = calInfo.correctTrialLength;                                                               %Correct Trial count
itCount = calInfo.incorrectTrialLength;                                                               %Incorrect Trial count
crCountC = calInfo.correctcorTrialLength;
crCountI = calInfo.incorrectcorTrialLength;
crcRow = [];                                                                %Correct correction trial row index
criRow = [];                                                                %Incorrect correction trial row index
cri = 1;                                                                    %Correct correction trial index marker
crc = 1;                                                                    %Incorrect correction trial index marker
dc = calInfo.DcorrectStartInd;
di = calInfo.DincorrectStartInd;
dcc = calInfo.DcorrectcorStartInd;
dic = calInfo.DincorrectcorStartInd;
dFrames = calInfo.DTrialLength;
dccount = 1;
dicount = 1;
dcccount = 1;
diccount = 1;
Dcoef = zeros(length(ms.FiltTraces(1,:)),1);
Acoef = Dcoef;

%Discecting the calcium traces in their respective trials
trace = transpose(normalize(ms.FiltTraces,'zscore'));
binarize = Binarize(ms);
btrace = transpose(binarize.binarizedTraces);

%Delay
dctrace = nan(length(ms.FiltTraces(1,:)),dFrames,length(dc));
ditrace = nan(length(ms.FiltTraces(1,:)),dFrames,length(di));
dcctrace = nan(length(ms.FiltTraces(1,:)),dFrames,length(dcc));

%Full traces
Ctrace = nan(length(ms.FiltTraces(1,:)),max(ctCount),length(ctInd(:,1)));   %Correct Trace
Itrace = nan(length(ms.FiltTraces(1,:)),max(itCount),length(itInd(:,1)));   %Incorrect Trace
CRCtrace = nan(length(ms.FiltTraces(1,:)),max(crCountC),length(crcInd(:,1)));   %Correct correction trace

if ~isempty(criRow)
    dicctrace = nan(length(ms.FiltTraces(1,:)),dFrames,length(dic));
    CRItrace = nan(length(ms.FiltTraces(1,:)),max(crCountI),length(criInd(:,1)));   %Incorrect correction trace
else
    CRItrace = [];
end

%extract respective traces from calcium
%Delay
dctrace = CutTraces(trace,dctrace,dc,dFrames*ones(length(dc),1),1);%Changed +60 to +dFrames
ditrace = CutTraces(trace,ditrace,di,dFrames*ones(length(di),1),1);
dcctrace = CutTraces(trace,dcctrace,dcc,dFrames*ones(length(dcc),1),1);

%Front Anchored
fCtrace = MatchLengths(CutTraces(trace,Ctrace,ctInd,ctCount,2),ctCount,2,[],[],[]);
fItrace = MatchLengths(CutTraces(trace,Itrace,itInd,itCount,2),itCount,2,[],[],[]);
fCRCtrace = MatchLengths(CutTraces(trace,CRCtrace,crcInd,crCountC,2),crCountC,2,[],[],[]);

%Back Anchored
bCtrace = MatchLengths(CutTraces(trace,Ctrace,ctInd,ctCount,3),ctCount,3,ctInd,dc,dFrames);
bItrace = MatchLengths(CutTraces(trace,Itrace,itInd,itCount,3),itCount,3,itInd,di,dFrames);
bCRCtrace = MatchLengths(CutTraces(trace,CRCtrace,crcInd,crCountC,3),crCountC,3,crcInd,dcc,dFrames);

if ~isempty(criRow)
    dicctrace = CutTraces(trace,dicctrace,dic,dFrames*ones(length(dic),1),1);
    fCRItrace = MatchLengths(CutTraces(trace,CRItrace,criInd,crCountI,2),crCountI,2,[],[],[]);
    bCRItrace = MatchLengths(CutTraces(trace,CRItrace,criInd,crCountI,3),crCountI,3,criInd,dic,dFrames);
else
    bCRItrace = [];
    fCRItrace = [];
    dicctrace = [];
end

%Shuffling
parfor i = 1: 100
    t = permute(trace, [2 1]);
    shuffledtraces = CShuffle(t);
    shuffledtraces = permute(shuffledtraces, [2 1]);
    icc = ~isempty(criRow);
    %Delay
    dshuffledC(:,i) = FullShuffleSplit(CutTraces(shuffledtraces,dctrace,dc,dFrames*ones(length(dc),1),1));%Changed +60 to +dFrames
    dshuffledI(:,i) = FullShuffleSplit(CutTraces(shuffledtraces,ditrace,di,dFrames*ones(length(di),1),1));
    dshuffledCcor(:,i) = FullShuffleSplit(CutTraces(shuffledtraces,dcctrace,dcc,dFrames*ones(length(dcc),1),1));
    
    %Front Anchored
    fshuffledC(:,i) = FullShuffleSplit(MatchLengths(CutTraces(shuffledtraces,Ctrace,ctInd,ctCount,2),ctCount,2,[],[],[]));
    fshuffledI(:,i) = FullShuffleSplit(MatchLengths(CutTraces(shuffledtraces,Itrace,itInd,itCount,2),itCount,2,[],[],[]));
    fshuffledCcor(:,i) = FullShuffleSplit(MatchLengths(CutTraces(shuffledtraces,CRCtrace,crcInd,crCountC,2),crCountC,2,[],[],[]));
    
    %Back Anchored
    bshuffledC(:,i) = FullShuffleSplit(MatchLengths(CutTraces(shuffledtraces,Ctrace,ctInd,ctCount,3),ctCount,3,ctInd,dc,dFrames));%Changed +60 to +dFrames
    bshuffledI(:,i) = FullShuffleSplit(MatchLengths(CutTraces(shuffledtraces,Itrace,itInd,itCount,3),itCount,3,itInd,di,dFrames));
    bshuffledCcor(:,i) = FullShuffleSplit(MatchLengths(CutTraces(shuffledtraces,CRCtrace,crcInd,crCountC,3),crCountC,3,crcInd,dcc,dFrames));
    
    %Incorrect Correction
    dtemp = zeros(ms.numNeurons,1);
    ftemp = zeros(ms.numNeurons,1);
    btemp = zeros(ms.numNeurons,1);
    if icc
        dtemp = FullShuffleSplit(CutTraces(shuffledtraces,dicctrace,dic,dFrames*ones(length(dic),1),1),dFrames,1);
        ftemp = FullShuffleSplit(MatchLengths(CutTraces(shuffledtraces,CRItrace,criInd,crCountI,2),crCountI,2,[],[],[]));
        btemp = FullShuffleSplit(MatchLengths(CutTraces(shuffledtraces,CRItrace,criInd,crCountI,3),crCountI,3,criInd,dic,dFrames));
    end
    
    dshuffledIcor(:,i) = dtemp;
    fshuffledIcor(:,i) = ftemp;
    bshuffledIcor(:,i) = btemp;
end

%Split-half reliability with random trial halves
out.DelaySplithalfcorrect = singleSplitShuffle(dctrace);
out.DelaySplithalfincorrect = singleSplitShuffle(ditrace);
out.DelaySplithalfccor = singleSplitShuffle(dcctrace);
out.DelaySplithalficor = singleSplitShuffle(dicctrace);

out.FrontSplithalfcorrect = singleSplitShuffle(fCtrace);
out.FrontSplithalfincorrect = singleSplitShuffle(fItrace);
out.FrontSplithalfccor = singleSplitShuffle(fCRCtrace);
out.FrontSplithalficor = singleSplitShuffle(fCRItrace);

out.BackSplithalfcorrect = singleSplitShuffle(bCtrace);
out.BackSplithalfincorrect = singleSplitShuffle(bItrace);
out.BackSplithalfccor = singleSplitShuffle(bCRCtrace);
out.BackSplithalficor = singleSplitShuffle(bCRItrace);

%Save shuffled correlation values
out.ScorrDelaycorrect = dshuffledC;
out.ScorrDelayincorrect = dshuffledI;
out.ScorrDelayccor = dshuffledCcor;
out.ScorrDelayicor = dshuffledIcor;

out.ScorrFrontcorrect= fshuffledC;
out.ScorrFrontincorrect = fshuffledI;
out.ScorrFrontccor = fshuffledCcor;
out.ScorrFronticor = fshuffledIcor;

out.ScorrBackcorrect = bshuffledC;
out.ScorrBackincorrect = bshuffledI;
out.ScorrBackccor = bshuffledCcor;
out.ScorrBackicor = bshuffledIcor;

%Find/Save cells which pass criteria
%Delay
out.Dpercentdc = prctile(out.ScorrDelaycorrect,99,2);                       % 99th percentile
out.DcorrectNumberPassed = length(find(out.DelaySplithalfcorrect> out.Dpercentdc));
out.DcorrectCells = find(out.DelaySplithalfcorrect> out.Dpercentdc);
out.DcorrectPassed = dctrace(out.DcorrectCells,:,:);
for i = length(out.DcorrectCells) : -1 : 1
    for j = 1 : length(calInfo.DcorrectStartInd)
        if any(find(btrace(out.DcorrectCells(i),calInfo.DcorrectStartInd(j):calInfo.DcorrectStartInd(j)+calInfo.DTrialLength)))
            break
        elseif j == length(calInfo.DcorrectStartInd)
            out.DcorrectCells(i) = [];
            out.DcorrectPassed(i,:,:) = [];
            out.DcorrectNumberPassed = out.DcorrectNumberPassed -1;
        end
    end
end

out.Dpercentdi = prctile(out.ScorrDelayincorrect,99,2);
out.DincorrectNumberPassed = length(find(out.DelaySplithalfincorrect> out.Dpercentdi));
out.DincorrectCells = find(out.DelaySplithalfincorrect> out.Dpercentdi);
out.DincorrectPassed = ditrace(out.DincorrectCells,:,:);
for i = length(out.DincorrectCells) : -1 : 1
    for j = 1 : length(calInfo.DincorrectStartInd)
        if any(find(btrace(out.DincorrectCells(i),calInfo.DincorrectStartInd(j):calInfo.DincorrectStartInd(j)+calInfo.DTrialLength)))
            break
        elseif j == length(calInfo.DincorrectStartInd)
            out.DincorrectCells(i) = [];
            out.DincorrectPassed(i,:,:) = [];
            out.DincorrectNumberPassed = out.DincorrectNumberPassed -1;
        end
    end
end

out.Dpercentdcc = prctile(out.ScorrDelayccor,99,2);
out.DccorNumberPassed = length(find(out.DelaySplithalfccor> out.Dpercentdcc));
out.DccorCells = find(out.DelaySplithalfccor> out.Dpercentdcc);
out.DccorPassed = dcctrace(out.DccorCells,:,:);
for i = length(out.DccorCells) : -1 : 1
    for j = 1 : length(calInfo.DcorrectcorStartInd)
        if any(find(btrace(out.DccorCells(i),calInfo.DcorrectcorStartInd(j):calInfo.DcorrectcorStartInd(j)+calInfo.DTrialLength)))
            break
        elseif j == length(calInfo.DcorrectcorStartInd)
            out.DccorCells(i) = [];
            out.DccorPassed(i,:,:) = [];
            out.DccorNumberPassed = out.DccorNumberPassed -1;
        end
    end
end

out.Dpercentdic = prctile(out.ScorrDelayicor,99,2);
out.DicorNumberPassed = length(find(out.DelaySplithalficor> out.Dpercentdic));
out.DicorCells = find(out.DelaySplithalficor> out.Dpercentdic);
out.DicorPassed = dicctrace(out.DicorCells,:,:);
for i = length(out.DicorCells) : -1 : 1
    for j = 1 : length(calInfo.DincorrectcorStartInd)
        if any(find(btrace(out.DincorrectcorCells(i),calInfo.DincorrectcorStartInd(j):calInfo.DincorrectcorStartInd(j)+calInfo.DTrialLength)))
            break
        elseif j == length(calInfo.DincorrectcorStartInd)
            out.DicorCells(i) = [];
            out.DicorPassed(i,:,:) = [];
            out.DicorNumberPassed = out.DicorNumberPassed -1;
        end
    end
end

%Front
out.Fpercentdc = prctile(out.ScorrFrontcorrect,99,2);
out.FcorrectNumberPassed = length(find(out.FrontSplithalfcorrect> out.Fpercentdc));
out.FcorrectCells = find(out.FrontSplithalfcorrect> out.Fpercentdc);
out.FcorrectPassed = fCtrace(out.FcorrectCells,:,:);
for i = length(out.FcorrectCells) : -1 : 1
    for j = 1 : length(calInfo.correctStartInd)
        if any(find(btrace(out.FcorrectCells(i),calInfo.correctStartInd(j):calInfo.correctStartInd(j)+calInfo.correctTrialLength(j))))
            break
        elseif j == length(calInfo.correctStartInd)
            out.FcorrectCells(i) = [];
            out.FcorrectPassed(i,:,:) = [];
            out.FcorrectNumberPassed = out.FcorrectNumberPassed -1;
        end
    end
end

out.Fpercentdi = prctile(out.ScorrFrontincorrect,99,2);
out.FincorrectNumberPassed = length(find(out.FrontSplithalfincorrect> out.Fpercentdi));
out.FincorrectCells = find(out.FrontSplithalfincorrect> out.Fpercentdi);
out.FincorrectPassed = fItrace(out.FincorrectCells,:,:);
for i = length(out.FincorrectCells) : -1 : 1
    for j = 1 : length(calInfo.incorrectStartInd)
        if any(find(btrace(out.FincorrectCells(i),calInfo.incorrectStartInd(j):calInfo.incorrectStartInd(j)+calInfo.incorrectTrialLength(j))))
            break
        elseif j == length(calInfo.incorrectStartInd)
            out.FincorrectCells(i) = [];
            out.FincorrectPassed(i,:,:) = [];
            out.FincorrectNumberPassed = out.FincorrectNumberPassed -1;
        end
    end
end

out.Fpercentdcc = prctile(out.ScorrFrontccor,99,2);
out.FccorNumberPassed = length(find(out.FrontSplithalfccor> out.Fpercentdcc));
out.FccorCells = find(out.FrontSplithalfccor> out.Fpercentdcc);
out.FccorPassed = fCRCtrace(out.FccorCells,:,:);
for i = length(out.FccorCells) : -1 : 1
    for j = 1 : length(calInfo.correctcorStartInd)
        if any(find(btrace(out.FccorCells(i),calInfo.correctcorStartInd(j):calInfo.correctcorStartInd(j)+calInfo.correctcorTrialLength(j))))
            break
        elseif j == length(calInfo.correctcorStartInd)
            out.FccorCells(i) = [];
            out.FccorPassed(i,:,:) = [];
            out.FccorNumberPassed = out.FccorNumberPassed -1;
        end
    end
end

out.Fpercentdic = prctile(out.ScorrFronticor,99,2);
out.FicorNumberPassed = length(find(out.FrontSplithalficor> out.Fpercentdic));
out.FicorCells = find(out.FrontSplithalficor> out.Fpercentdic);
out.FicorPassed = fCRItrace(out.FicorCells,:,:);
for i = length(out.FicorCells) : -1 : 1
    for j = 1 : length(calInfo.incorrectcorStartInd)
        if any(find(btrace(out.FicorCells(i),calInfo.incorrectcorStartInd(j):calInfo.incorrectcorStartInd(j)+calInfo.incorrectcorTrialLength(j))))
            break
        elseif j == length(calInfo.incorrectcorStartInd)
            out.FicorCells(i) = [];
            out.FicorPassed(i,:,:) = [];
            out.FicorNumberPassed = out.FicorNumberPassed -1;
        end
    end
end

%Back
out.Bpercentdc = prctile(out.ScorrBackcorrect,99,2);
out.BcorrectNumberPassed = length(find(out.BackSplithalfcorrect> out.Bpercentdc));
out.BcorrectCells = find(out.BackSplithalfcorrect> out.Bpercentdc);
out.BcorrectPassed = bCtrace(out.BcorrectCells,:,:);
for i = length(out.BcorrectCells) : -1 : 1
    for j = 1 : length(calInfo.correctStartInd)
        if any(find(btrace(out.BcorrectCells(i),calInfo.correctStartInd(j):calInfo.correctStartInd(j)+calInfo.correctTrialLength(j))))
            break
        elseif j == length(calInfo.correctStartInd)
            out.BcorrectCells(i) = [];
            out.BcorrectPassed(i,:,:) = [];
            out.BcorrectNumberPassed = out.BcorrectNumberPassed -1;
        end
    end
end

out.Bpercentdi = prctile(out.ScorrBackincorrect,99,2);
out.BincorrectNumberPassed = length(find(out.BackSplithalfincorrect> out.Bpercentdi));
out.BincorrectCells = find(out.BackSplithalfincorrect> out.Bpercentdi);
out.BincorrectPassed = bItrace(out.BincorrectCells,:,:);
for i = length(out.BincorrectCells) : -1 : 1
    for j = 1 : length(calInfo.incorrectStartInd)
        if any(find(btrace(out.BincorrectCells(i),calInfo.incorrectStartInd(j):calInfo.incorrectStartInd(j)+calInfo.incorrectTrialLength(j))))
            break
        elseif j == length(calInfo.incorrectStartInd)
            out.BincorrectCells(i) = [];
            out.BincorrectPassed(i,:,:) = [];
            out.BincorrectNumberPassed = out.BincorrectNumberPassed -1;
        end
    end
end

out.Bpercentdcc = prctile(out.ScorrBackccor,99,2);
out.BccorNumberPassed = length(find(out.BackSplithalfccor> out.Bpercentdcc));
out.BccorCells = find(out.BackSplithalfccor> out.Bpercentdcc);
out.BccorPassed = bCRCtrace(out.BccorCells,:,:);
for i = length(out.BccorCells) : -1 : 1
    for j = 1 : length(calInfo.correctcorStartInd)
        if any(find(btrace(out.BccorCells(i),calInfo.correctcorStartInd(j):calInfo.correctcorStartInd(j)+calInfo.correctcorTrialLength(j))))
            break
        elseif j == length(calInfo.correctcorStartInd)
            out.BccorCells(i) = [];
            out.BccorPassed(i,:,:) = [];
            out.BccorNumberPassed = out.BccorNumberPassed -1;
        end
    end
end

out.Bpercentdic = prctile(out.ScorrBackicor,99,2);
out.BicorNumberPassed = length(find(out.BackSplithalficor> out.Bpercentdic));
out.BicorCells = find(out.BackSplithalficor> out.Bpercentdic);
out.BicorPassed = bCRItrace(out.BicorCells,:,:);
for i = length(out.BicorCells) : -1 : 1
    for j = 1 : length(calInfo.incorrectcorStartInd)
        if any(find(btrace(out.BicorCells(i),calInfo.incorrectcorStartInd(j):calInfo.incorrectcorStartInd(j)+calInfo.incorrectcorTrialLength(j))))
            break
        elseif j == length(calInfo.incorrectcorStartInd)
            out.BicorCells(i) = [];
            out.BicorPassed(i,:,:) = [];
            out.BicorNumberPassed = out.BicorNumberPassed -1;
        end
    end    
end

toc
end

function rast = sortPeaks(rast)
[~,maxind] = max(rast,[],2);
[~, rastsort] = sort(maxind);
rast = rast(rastsort,:);
end

%This takes the split half reliabitily of the fully shuffled data
function [out] = FullShuffleSplit(trial)
if ~isempty(trial)
    if ~(length(trial(1,1,:)) == 0)
        out = zeros(length(trial(:,1,1)),1);
        if ~isempty(trial) && length(trial(1,1,:))>1
            temp = permute(trial,[3 2 1]);
%             valsum = sum(temp,2);%%%%%%%%%%%%%%%%%%%%%%%
            for i = 1 : length(temp(1,1,:))
                damp = 1;%length(find(valsum(:,1,i)))/length(temp(:,1,1));%%%%%%%%%%%%%%%%%
                h1 = nanmean(temp((1:round(length(temp(:,1))/2)),:,i),1);
                h2 = nanmean(temp((round(length(temp(:,1))/2)+1:end),:,i),1);
                out(i) = corr2(h1,h2)*damp;%%%%%%%%%%%%%%%%%%%
            end
        end
    else
        out = NaN;
    end
else
    out = NaN;
end
end

%Given the indices of each trial, will return a matrix of calcium within
%each trial
function [out] = CutTraces(traces,trial,ind, count,cut)
if cut == 2 || cut == 1
    for i =1: length(trial(1,1,:))
        trial(:,1:count(i),i) = traces(:,ind(i):ind(i)+count(i)-1);
    end
elseif cut == 3
    for i =1: length(trial(1,1,:))
        trial(:, (end- count(i))+1:end,i) = traces(:,ind(i):ind(i)+count(i)-1);
    end
end
out = trial;
end

% % % %Will match the length of all trials with the shortest trial length,
% % % %dropping any extra components
% % % function [out] = MatchLengths(trial,cl,cut,ind,dind,dframe)
% % % if ~isempty(trial)
% % %     if ~(length(trial(1,1,:)) == 0)
% % %         out = zeros(size(trial ));
% % %         if ~isempty(trial) && length(trial(1,1,:))>1
% % %             temp = permute(trial,[3 2 1]);
% % %             if cut == 2 %Front Anchored: eliminate back end excess
% % %                 if ~(min(temp) == length(temp))
% % %                     out = permute(temp(:,1:min(cl),:),[3 2 1]);
% % %                 end
% % %             elseif cut == 3 && ~isempty(ind) && ~isempty(dind) && ~isempty(dframe) %Back Anchored: eliminate front end excess
% % %                 sT = cl - ((dind(:) - ind(:))+dframe);                      %Shortest Trial for the BA is found after subtracting the delay period
% % %                 if ~(min(temp) == length(temp))
% % %                     out = permute(temp(:,end-min(sT)+1:end,:),[3 2 1]);
% % %                 end
% % %             end            
% % %         end
% % %     else
% % %         out = NaN;
% % %     end
% % % else
% % %     out = NaN;
% % % end    
% % % end