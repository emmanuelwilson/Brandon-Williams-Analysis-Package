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

function [out] = TouchScreenShuffleSeperation_V2Weighted(SepDelay,SepFront,SepBack,ms,calInfo,events)
%seperate frame map from events file
tic
for i = 1 : length(events)
    frameMap(i,1) = str2num(char(events(i,1)));
end
NumberOfChoices = 5;
%Initilize variables
out.ScorrDelaycorrect = zeros(length(ms.FiltTraces(1,:)),100,NumberOfChoices);
out.ScorrDelayincorrect = zeros(length(ms.FiltTraces(1,:)),100,NumberOfChoices);
out.ScorrDelayccor = zeros(length(ms.FiltTraces(1,:)),100,NumberOfChoices);
out.ScorrDelayicor = zeros(length(ms.FiltTraces(1,:)),100,NumberOfChoices);

out.DelaySplithalfcorrect = zeros(length(ms.FiltTraces(1,:)),NumberOfChoices);
out.DelaySplithalfincorrect = zeros(length(ms.FiltTraces(1,:)),NumberOfChoices);
out.DelaySplithalfccor = zeros(length(ms.FiltTraces(1,:)),NumberOfChoices);
out.DelaySplithalficor = zeros(length(ms.FiltTraces(1,:)),NumberOfChoices);

out.ScorrFrontcorrect = zeros(length(ms.FiltTraces(1,:)),100,NumberOfChoices);
out.ScorrFrontincorrect = zeros(length(ms.FiltTraces(1,:)),100,NumberOfChoices);
out.ScorrFrontccor = zeros(length(ms.FiltTraces(1,:)),100,NumberOfChoices);
out.ScorrFronticor = zeros(length(ms.FiltTraces(1,:)),100,NumberOfChoices);

out.FrontSplithalfcorrect = zeros(length(ms.FiltTraces(1,:)),NumberOfChoices);
out.FrontSplithalfincorrect = zeros(length(ms.FiltTraces(1,:)),NumberOfChoices);
out.FrontSplithalfccor = zeros(length(ms.FiltTraces(1,:)),NumberOfChoices);
out.FrontSplithalficor = zeros(length(ms.FiltTraces(1,:)),NumberOfChoices);

out.ScorrBackcorrect = zeros(length(ms.FiltTraces(1,:)),100,NumberOfChoices);
out.ScorrBackincorrect = zeros(length(ms.FiltTraces(1,:)),100,NumberOfChoices);
out.ScorrBackccor = zeros(length(ms.FiltTraces(1,:)),100,NumberOfChoices);
out.ScorrBackicor = zeros(length(ms.FiltTraces(1,:)),100,NumberOfChoices);

out.BackSplithalfcorrect = zeros(length(ms.FiltTraces(1,:)),NumberOfChoices);
out.BackSplithalfincorrect = zeros(length(ms.FiltTraces(1,:)),NumberOfChoices);
out.BackSplithalfccor = zeros(length(ms.FiltTraces(1,:)),NumberOfChoices);
out.BackSplithalficor = zeros(length(ms.FiltTraces(1,:)),NumberOfChoices);

ctInd = calInfo.correctStartInd;                                                                 %Correct Trial Indices
itInd = calInfo.incorrectStartInd;                                                                 %Incorrect Trial Indices
ctCount = calInfo.correctTrialLength;                                                               %Correct Trial count
itCount = calInfo.incorrectTrialLength;                                                               %Incorrect Trial count
crCountC = calInfo.correctcorTrialLength;
crCountI = calInfo.incorrectcorTrialLength;
dFrames = calInfo.DTrialLength;
crcRow = [];                                                                %Correct correction trial row index
criRow = [];                                                                %Incorrect correction trial row index
crcInd = calInfo.correctcorStartInd;                                                                %Correct correction trial indices
criInd = calInfo.incorrectcorStartInd;                                                                %Incorrect correction trial indices
cri = 1;                                                                    %Correct correction trial index marker
crc = 1;                                                                    %Incorrect correction trial index marker
dc = calInfo.DcorrectStartInd;
di = calInfo.DincorrectStartInd;
dcc = calInfo.DcorrectcorStartInd;
dic = calInfo.DincorrectcorStartInd;
dccount = 1;
dicount = 1;
dcccount = 1;
diccount = 1;

% trace = transpose(normalize(ms.FiltTraces,'zscore'));
trace = ms.FiltTraces';
binarize = Binarize(ms);
btrace = transpose(binarize.binarizedTraces);

%Shuffling
for i = 1: 100
    t = permute(trace, [2 1]);
    shuffledtraces = CShuffle(t);
    shuffledtraces = permute(shuffledtraces, [2 1]);
    for c = 1: NumberOfChoices        
%         for j = length(ms.FiltTraces(1,:))
            %Correct Trials
            if ~isempty(calInfo.Ctrial{c})                
                dshuffledC{c}(:,i) = FullShuffleSplit(CutTraces(shuffledtraces,permute(SepDelay.CorrectTrial{c},[3 2 1]),dc(calInfo.Ctrial{c}),dFrames*ones(length(dc(calInfo.Ctrial{c})),1),1));%Changed +60 to +dFrames                
                fshuffledC{c}(:,i) = FullShuffleSplit(CutTraces(shuffledtraces,permute(SepFront.CorrectTrial{c},[3 2 1]),ctInd(calInfo.Ctrial{c}),ctCount(calInfo.Ctrial{c}),2));
                bshuffledC{c}(:,i) = FullShuffleSplit(CutTraces(shuffledtraces,permute(SepBack.CorrectTrial{c},[3 2 1]),ctInd(calInfo.Ctrial{c}),ctCount(calInfo.Ctrial{c}),3));%Changed +60 to +dFrames
            else
                dshuffledC{c} = [];
                fshuffledC{c} = [];                
                bshuffledC{c} = [];
            end
            
            %Incorrect Trials
            if ~isempty(calInfo.Itrial{c})                
                dshuffledI{c}(:,i) = FullShuffleSplit(CutTraces(shuffledtraces,permute(SepDelay.IncorrectTrial{c},[3 2 1]),di(calInfo.Itrial{c}),dFrames*ones(length(di(calInfo.Itrial{c})),1),1));                
                fshuffledI{c}(:,i) = FullShuffleSplit(CutTraces(shuffledtraces,permute(SepFront.IncorrectTrial{c},[3 2 1]),itInd(calInfo.Itrial{c}),itCount(calInfo.Itrial{c}),2));
                bshuffledI{c}(:,i) = FullShuffleSplit(CutTraces(shuffledtraces,permute(SepBack.IncorrectTrial{c},[3 2 1]),itInd(calInfo.Itrial{c}),itCount(calInfo.Itrial{c}),3));
            else
                dshuffledI{c} = [];                
                fshuffledI{c} = [];                
                bshuffledI{c} = [];
            end
            
            %Correct Correction
            if ~isempty(calInfo.CCtrial{c})         
                dshuffledCcor{c}(:,i) = FullShuffleSplit(CutTraces(shuffledtraces,permute(SepDelay.CorrectCorrectionTrial{c},[3 2 1]),dcc(calInfo.CCtrial{c}),dFrames*ones(length(dcc(calInfo.CCtrial{c})),1),1));                                
                bshuffledCcor{c}(:,i) = FullShuffleSplit(CutTraces(shuffledtraces,permute(SepFront.CorrectCorrectionTrial{c},[3 2 1]),crcInd(calInfo.CCtrial{c}),crCountC(calInfo.CCtrial{c}),3));
                fshuffledCcor{c}(:,i) = FullShuffleSplit(CutTraces(shuffledtraces,permute(SepBack.CorrectCorrectionTrial{c},[3 2 1]),crcInd(calInfo.CCtrial{c}),crCountC(calInfo.CCtrial{c}),2));
            else
                dshuffledCcor{c} = [];
                fshuffledCcor{c} = [];
                bshuffledCcor{c} = [];
            end
            %Incorrect Correction
            if ~isempty(calInfo.ICtrial{c})
                dshuffledIcor{c}(:,i) = FullShuffleSplit(CutTraces(shuffledtraces,permute(SepDelay.IncorrectCorrectionTrial{c},[3 2 1]),dic(calInfo.ICtrial{c}),dFrames*ones(length(dic(calInfo.ICtrial{c})),1),1));
                fshuffledIcor{c}(:,i) = FullShuffleSplit(CutTraces(shuffledtraces,permute(SepFront.IncorrectCorrectionTrial{c},[3 2 1]),criInd(calInfo.ICtrial{c}),crCountI(calInfo.ICtrial{c}),2));
                bshuffledIcor{c}(:,i) = FullShuffleSplit(CutTraces(shuffledtraces,permute(SepBack.IncorrectCorrectionTrial{c},[3 2 1]),criInd(calInfo.ICtrial{c}),crCountI(calInfo.ICtrial{c}),3));
            else
                fshuffledIcor{c} = [];
                dshuffledIcor{c} = [];
                bshuffledIcor{c} = [];
            end
%         end
    end
end

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

%Split-half reliability with random trial halves
for c = 1 : NumberOfChoices
    %Correct
    if ~isempty(calInfo.Ctrial{c})
        out.DelaySplithalfcorrect(:,c) = singleSplitShuffleWeighted(permute(SepDelay.CorrectTrial{c},[3 2 1]));
        out.FrontSplithalfcorrect(:,c) = singleSplitShuffleWeighted(permute(SepFront.CorrectTrial{c},[3 2 1]));
        out.BackSplithalfcorrect(:,c) = singleSplitShuffleWeighted(permute(SepBack.CorrectTrial{c},[3 2 1]));
        out.Dpercentdc(:,c) = prctile(out.ScorrDelaycorrect{c},99,2);                       % 99th percentile
        
        %Delay
        out.DcorrectNumberPassed{c} = length(find(out.DelaySplithalfcorrect(:,c)> out.Dpercentdc(:,c)));
        out.DcorrectCells{c} = find(out.DelaySplithalfcorrect(:,c)> out.Dpercentdc(:,c));
        out.DcorrectPassed{c} = SepDelay.CorrectTrial{c}(:,:,out.DcorrectCells{c});
        for i = length(out.DcorrectCells{c}) : -1 : 1
            for j = 1 : length(calInfo.Ctrial{c})
                if any(find(btrace(out.DcorrectCells{c}(i),calInfo.DcorrectStartInd(calInfo.Ctrial{c}(j)):calInfo.DcorrectStartInd(calInfo.Ctrial{c}(j))+calInfo.DTrialLength)))
                    break
                elseif j == length(calInfo.Ctrial{c})
                    out.DcorrectCells{c}(i) = [];
                    out.DcorrectPassed{c}(:,:,i) = [];
                    out.DcorrectNumberPassed{c} = out.DcorrectNumberPassed{c} -1;
                end
            end
        end
        
        %Front
        out.Fpercentdc(:,c) = prctile(out.ScorrFrontcorrect{c},99,2);
        out.FcorrectNumberPassed{c} = length(find(out.FrontSplithalfcorrect(:,c)> out.Fpercentdc(:,c)));
        out.FcorrectCells{c} = find(out.FrontSplithalfcorrect(:,c)> out.Fpercentdc(:,c));
        out.FcorrectPassed{c} = SepFront.CorrectTrial{c}(:,:,out.FcorrectCells{c});
        for i = length(out.FcorrectCells{c}) : -1 : 1
            for j = 1 : length(calInfo.Ctrial{c})
                if any(find(btrace(out.FcorrectCells{c}(i),calInfo.correctStartInd(calInfo.Ctrial{c}(j)):calInfo.correctStartInd(calInfo.Ctrial{c}(j))+calInfo.correctTrialLength(calInfo.Ctrial{c}(j)))))
                    break
                elseif j == length(calInfo.Ctrial{c})
                    out.FcorrectCells{c}(i) = [];
                    out.FcorrectPassed{c}(:,:,i) = [];
                    out.FcorrectNumberPassed{c} = out.FcorrectNumberPassed{c} -1;
                end
            end
        end
        
        %Back
        out.Bpercentdc(:,c) = prctile(out.ScorrBackcorrect{c},99,2);
        out.BcorrectNumberPassed{c} = length(find(out.BackSplithalfcorrect(:,c)> out.Bpercentdc(:,c)));
        out.BcorrectCells{c} = find(out.BackSplithalfcorrect(:,c)> out.Bpercentdc(:,c));
        out.BcorrectPassed{c} = SepBack.CorrectTrial{c}(:,:,out.BcorrectCells{c});
        for i = length(out.BcorrectCells{c}) : -1 : 1
            for j = 1 : length(calInfo.Ctrial{c})
                if any(find(btrace(out.BcorrectCells{c}(i),calInfo.correctStartInd(calInfo.Ctrial{c}(j)):calInfo.correctStartInd(calInfo.Ctrial{c}(j))+calInfo.correctTrialLength(calInfo.Ctrial{c}(j)))))
                    break
                elseif j == length(calInfo.Ctrial{c})
                    out.BcorrectCells{c}(i) = [];
                    out.BcorrectPassed{c}(:,:,i) = [];
                    out.BcorrectNumberPassed{c} = out.BcorrectNumberPassed{c} -1;
                end
            end
        end
    end
    
    
    %Incorrect
    if ~isempty(calInfo.Itrial{c})
        out.DelaySplithalfincorrect(:,c) = singleSplitShuffleWeighted(permute(SepDelay.IncorrectTrial{c},[3 2 1]));
        out.FrontSplithalfincorrect(:,c) = singleSplitShuffleWeighted(permute(SepFront.IncorrectTrial{c},[3 2 1]));
        out.BackSplithalfincorrect(:,c) = singleSplitShuffleWeighted(permute(SepBack.IncorrectTrial{c},[3 2 1]));
        
        %Delay
        out.Dpercentdi(:,c) = prctile(out.ScorrDelayincorrect{c},99,2);                       % 99th percentile
        out.DincorrectNumberPassed{c} = length(find(out.DelaySplithalfincorrect(:,c)> out.Dpercentdi(:,c)));
        out.DincorrectCells{c} = find(out.DelaySplithalfincorrect(:,c)> out.Dpercentdi(:,c));
        out.DincorrectPassed{c} = SepDelay.IncorrectTrial{c}(:,:,out.DincorrectCells{c});
        for i = length(out.DincorrectCells{c}) : -1 : 1            
            for j = 1 : length(calInfo.Itrial{c})
                if any(find(btrace(out.DincorrectCells{c}(i),calInfo.DincorrectStartInd(calInfo.Itrial{c}(j)):calInfo.DincorrectStartInd(calInfo.Itrial{c}(j))+calInfo.DTrialLength)))
                    break
                elseif j == length(calInfo.Itrial{c})
                    out.DincorrectCells{c}(i) = [];
                    out.DincorrectPassed{c}(:,:,i) = [];
                    out.DincorrectNumberPassed{c} = out.DincorrectNumberPassed{c} -1;
                end
            end
        end
        
        %Front
        out.Fpercentdi(:,c) = prctile(out.ScorrFrontincorrect{c},99,2);
        out.FincorrectNumberPassed{c} = length(find(out.FrontSplithalfincorrect(:,c)> out.Fpercentdi(:,c)));
        out.FincorrectCells{c} = find(out.FrontSplithalfincorrect(:,c)> out.Fpercentdi(:,c));
        out.FincorrectPassed{c} = SepFront.IncorrectTrial{c}(:,:,out.FincorrectCells{c});
        for i = length(out.FincorrectCells{c}) : -1 : 1
            for j = 1 : length(calInfo.Itrial{c})
                if any(find(btrace(out.FincorrectCells{c}(i),calInfo.incorrectStartInd(calInfo.Itrial{c}(j)):calInfo.incorrectStartInd(calInfo.Itrial{c}(j))+calInfo.incorrectTrialLength(calInfo.Itrial{c}(j)))))
                    break
                elseif j == length(calInfo.Itrial{c})
                    out.FincorrectCells{c}(i) = [];
                    out.FincorrectPassed{c}(:,:,i) = [];
                    out.FincorrectNumberPassed{c} = out.FincorrectNumberPassed{c} -1;
                end
            end
        end
        
        %Back
        out.Bpercentdi(:,c) = prctile(out.ScorrBackincorrect{c},99,2);
        out.BincorrectNumberPassed{c} = length(find(out.BackSplithalfincorrect(:,c)> out.Bpercentdi(:,c)));
        out.BincorrectCells{c} = find(out.BackSplithalfincorrect(:,c)> out.Bpercentdi(:,c));
        out.BincorrectPassed{c} = SepBack.IncorrectTrial{c}(:,:,out.BincorrectCells{c});
        for i = length(out.BincorrectCells{c}) : -1 : 1
            for j = 1 : length(calInfo.Itrial{c})
                if any(find(btrace(out.BincorrectCells{c}(i),calInfo.incorrectStartInd(calInfo.Itrial{c}(j)):calInfo.incorrectStartInd(calInfo.Itrial{c}(j))+calInfo.incorrectTrialLength(calInfo.Itrial{c}(j)))))
                    break
                elseif j == length(calInfo.Itrial{c})
                    out.BincorrectCells{c}(i) = [];
                    out.BincorrectPassed{c}(:,:,i) = [];
                    out.BincorrectNumberPassed{c} = out.BincorrectNumberPassed{c} -1;
                end
            end
        end
    end
    
    %Correct Correction
    if ~isempty(calInfo.CCtrial{c})
        out.DelaySplithalfccor(:,c) = singleSplitShuffleWeighted(permute(SepDelay.CorrectCorrectionTrial{c},[3 2 1]));
        out.FrontSplithalfccor(:,c) = singleSplitShuffleWeighted(permute(SepFront.CorrectCorrectionTrial{c},[3 2 1]));
        out.BackSplithalfccor(:,c) = singleSplitShuffleWeighted(permute(SepBack.CorrectCorrectionTrial{c},[3 2 1]));
        
        %Delay
        out.Dpercentdcc(:,c) = prctile(out.ScorrDelayccor{c},99,2);
        out.DccorNumberPassed{c} = length(find(out.DelaySplithalfccor(:,c)> out.Dpercentdcc(:,c)));
        out.DccorCells{c} = find(out.DelaySplithalfccor(:,c)> out.Dpercentdcc(:,c));
        out.DccorPassed{c} = SepDelay.CorrectCorrectionTrial{c}(:,:,out.DccorCells{c});
        for i = length(out.DccorCells{c}) : -1 : 1
            for j = 1 : length(calInfo.CCtrial{c})
                if any(find(btrace(out.DccorCells{c}(i),calInfo.DcorrectcorStartInd(calInfo.CCtrial{c}(j)):calInfo.DcorrectcorStartInd(calInfo.CCtrial{c}(j))+calInfo.DTrialLength)))
                    break
                elseif j == length(calInfo.CCtrial{c})
                    out.DccorCells{c}(i) = [];
                    out.DccorPassed{c}(:,:,i) = [];
                    out.DccorNumberPassed{c} = out.DccorNumberPassed{c} -1;
                end
            end
        end
        
        %Front
        out.Fpercentdcc(:,c) = prctile(out.ScorrFrontccor{c},99,2);
        out.FccorNumberPassed{c} = length(find(out.FrontSplithalfccor(:,c)> out.Fpercentdcc(:,c)));
        out.FccorCells{c} = find(out.FrontSplithalfccor(:,c)> out.Fpercentdcc(:,c));
        out.FccorPassed{c} = SepFront.CorrectCorrectionTrial{c}(:,:,out.FccorCells{c});
        for i = length(out.FccorCells{c}) : -1 : 1
            if ~isempty(calInfo.CCtrial{c})
                for j = 1 : length(calInfo.CCtrial{c})
                    if any(find(btrace(out.FccorCells{c}(i),calInfo.correctcorStartInd(calInfo.CCtrial{c}(j)):calInfo.correctcorStartInd(calInfo.CCtrial{c}(j))+calInfo.correctcorTrialLength(calInfo.CCtrial{c}(j)))))
                        break
                    elseif j == length(calInfo.CCtrial{c})
                        out.FccorCells{c}(i) = [];
                        out.FccorPassed{c}(:,:,i) = [];
                        out.FccorNumberPassed{c} = out.FccorNumberPassed{c} -1;
                    end
                end
            end
        end
        
        %Back
        out.Bpercentdcc(:,c) = prctile(out.ScorrBackccor{c},99,2);
        out.BccorNumberPassed{c} = length(find(out.BackSplithalfccor(:,c)> out.Bpercentdcc(:,c)));
        out.BccorCells{c} = find(out.BackSplithalfccor(:,c)> out.Bpercentdcc(:,c));
        out.BccorPassed{c} = SepBack.CorrectCorrectionTrial{c}(:,:,out.BccorCells{c});
        for i = length(out.BccorCells{c}) : -1 : 1
            for j = 1 : length(calInfo.CCtrial{c})
                if any(find(btrace(out.BccorCells{c}(i),calInfo.correctcorStartInd(calInfo.CCtrial{c}(j)):calInfo.correctcorStartInd(calInfo.CCtrial{c}(j))+calInfo.correctcorTrialLength(calInfo.CCtrial{c}(j)))))
                    break
                elseif j == length(calInfo.CCtrial{c})
                    out.BccorCells{c}(i) = [];
                    out.BccorPassed{c}(:,:,i) = [];
                    out.BccorNumberPassed{c} = out.BccorNumberPassed{c} -1;
                end
            end
        end
    end
    %Incorrect Correction
    if ~isempty(calInfo.ICtrial{c})
        out.DelaySplithalficor(:,c) = singleSplitShuffleWeighted(permute(SepDelay.IncorrectCorrectionTrial{c},[3 2 1]));
        out.FrontSplithalficor(:,c) = singleSplitShuffleWeighted(permute(SepFront.IncorrectCorrectionTrial{c},[3 2 1]));
        out.BackSplithalficor(:,c) = singleSplitShuffleWeighted(permute(SepBack.IncorrectCorrectionTrial{c},[3 2 1]));
        
        %Delay
        out.Dpercentdic(:,c) = prctile(out.ScorrDelayicor{c},99,2);
        out.DicorNumberPassed{c} = length(find(out.DelaySplithalficor(:,c)> out.Dpercentdic(:,c)));
        out.DicorCells{c} = find(out.DelaySplithalficor(:,c)> out.Dpercentdic(:,c));
        out.DicorPassed{c} = SepDelay.IncorrectCorrectionTrial{c}(:,:,out.DicorCells{c});
        for i = length(out.DicorCells{c}) : -1 : 1
            for j = 1 : length(calInfo.ICtrial{c})
                if any(find(btrace(out.DicorCells{c}(i),calInfo.DincorrectcorStartInd(calInfo.ICtrial{c}(j)):calInfo.DincorrectcorStartInd(calInfo.ICtrial{c}(j))+calInfo.DTrialLength)))
                    break
                elseif j == length(calInfo.ICtrial{c})
                    out.DicorCells{c}(i) = [];
                    out.DicorPassed{c}(:,:,i) = [];
                    out.DicorNumberPassed{c} = out.DicorNumberPassed{c} -1;
                end
            end
        end
        
        %Front
        out.Fpercentdic(:,c) = prctile(out.ScorrFronticor{c},99,2);
        out.FicorNumberPassed{c} = length(find(out.FrontSplithalficor(:,c)> out.Fpercentdic(:,c)));
        out.FicorCells{c} = find(out.FrontSplithalficor(:,c)> out.Fpercentdic(:,c));
        out.FicorPassed{c} = SepFront.IncorrectCorrectionTrial{c}(:,:,out.FicorCells{c});
        for i = length(out.FicorCells{c}) : -1 : 1
            for j = 1 : length(calInfo.ICtrial{c})
                if any(find(btrace(out.FicorCells{c}(i),calInfo.incorrectcorStartInd(calInfo.ICtrial{c}(j)):calInfo.incorrectcorStartInd(calInfo.ICtrial{c}(j))+calInfo.incorrectcorTrialLength(calInfo.ICtrial{c}(j)))))
                    break
                elseif j == length(calInfo.ICtrial{c})
                    out.FicorCells{c}(i) = [];
                    out.FicorPassed{c}(:,:,i) = [];
                    out.FicorNumberPassed{c} = out.FicorNumberPassed{c} -1;
                end
            end
        end
        
        %Back
        out.Bpercentdic(:,c) = prctile(out.ScorrBackicor{c},99,2);
        out.BicorNumberPassed{c} = length(find(out.BackSplithalficor(:,c)> out.Bpercentdic(:,c)));
        out.BicorCells{c} = find(out.BackSplithalficor(:,c)> out.Bpercentdic(:,c));
        out.BicorPassed{c} = SepBack.IncorrectCorrectionTrial{c}(:,:,out.BicorCells{c});
        for i = length(out.BicorCells{c}) : -1 : 1
            for j = 1 : length(calInfo.ICtrial{c})
                if any(find(btrace(out.BicorCells{c}(i),calInfo.incorrectcorStartInd(calInfo.ICtrial{c}(j)):calInfo.incorrectcorStartInd(calInfo.ICtrial{c}(j))+calInfo.incorrectcorTrialLength(calInfo.ICtrial{c}(j)))))
                    break
                elseif j == length(calInfo.ICtrial{c})
                    out.BicorCells{c}(i) = [];
                    out.BicorPassed{c}(:,:,i) = [];
                    out.BicorNumberPassed{c} = out.BicorNumberPassed{c} -1;
                end
            end
            
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
            valsum = sum(temp,2);               
            count = 1;
            for i = 1 : length(temp(1,1,:))
                damp = length(find(valsum(:,1,i)))/length(temp(:,1,1));
                h1 = nanmean(temp((1:round(length(temp(:,1))/2)),:,i),1);
                h2 = nanmean(temp((round(length(temp(:,1))/2)+1:end),:,i),1);
                out(i) = corr2(h1,h2)*damp;
                if out(i) > 0.9
                end
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
if cut == 1
    for i =1: length(trial(1,1,:))
        trial(:,:,i) = traces(:,ind(i):ind(i)+min(count));
    end
elseif cut == 3 || cut == 2
    for i =1: length(trial(1,1,:))
        trial(:,:,i) = traces(:,ind(i):ind(i)+length(trial(1,:,1))-1);
    end
end
out = trial;
end