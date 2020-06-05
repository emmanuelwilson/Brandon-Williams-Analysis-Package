%% Find Confidence interval of Mean Fluorescence activity
%
%INPUT: -rasts: Structure which contains calcium responses of cells of all
%           trials. Each trial type is saved in a cell array with each
%           column corresponding to the choice location made by the mouse. 
%       -cut:
%           -1: Delay anchor
%           -2: FrontAnchor
%           -3: BackAnchor
%       -path: Folder path of current session being analyzed (where you
%           want to save your results)
%OUTPUT:Saves all figures in respective folder (Delay, Front, etc) and
%           matrices of said figures. 
%       -MatStruct: Contains matrices to replicate saved figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Emmanuel Wilson

function out = ErrorSource_V2_invisible(rasts,cut,path,indpassed)
fps = 30;
ticmultiplier = 2;
%Find/create right folder path to save results
if cut == 1
    path = [path,'/Delay'];
elseif cut == 2
    path = [path,'/FrontAnchored'];
elseif cut == 3
    path = [path,'/BackAnchored'];
end
%Choice Counters
m1 = 0;
m2 = 0;
%Find number of choices
for i = 1: length(rasts.CorrectTrial)
    m1 = m1 + ~isempty(rasts.CorrectTrial{i});
    m1i(i) = ~isempty(rasts.CorrectTrial{i});
    m2 = m2 + ~isempty(rasts.IncorrectTrial{i});
    m2i(i) = ~isempty(rasts.IncorrectTrial{i});
end

%Max number of choices made
maxval = max(cat(1,m1,m2)); 

%Find which trial had the most choices
if m1 == maxval
    poke = find(m1i);
elseif m2 == maxval
    poke = find(m2i);
end
f = figure(1);
f.Visible = 'off';
%Plot Generation and data processing
for i = 1: maxval
    pos = poke(i);                                                         %Find position of choice
    tempC = ExtractTrial(rasts,1,pos);                                     %Extract Correct traces
    tempI = ExtractTrial(rasts,2,pos);                                     %Extract Incorrect traces
    mkdir(path,['ConfidenceIntervals/',num2str(pos)])                      %Make folder to save figures/results
    
    %find start/end frames for comparison loop
    if ~isempty(tempI) && ~isempty(tempC) && length(tempC(1,:,1)) > length(tempI(1,:,1))
        if cut == 1 || cut == 2
            timeS = 1;
            timeE = length(tempI(1,:,1));
        elseif cut == 3
            timeS = length(tempC(1,:,1))-length(tempI(1,:,1));
            timeE = length(tempC(1,:,1));
            addon = NaN(length(tempI(:,1,1)),timeS,length(tempI(1,1,:)));
            tempI = cat(2,addon,tempI);
        end
    elseif (isempty(tempI) && ~isempty(tempC)) || length(tempI(1,:,1)) > length(tempC(1,:,1))
        if cut == 1 || cut == 2
            timeS = 1;
            timeE = length(tempC(1,:,1));
        elseif cut == 3
            if (isempty(tempI) && ~isempty(tempC)) 
                timeS = 1;
                timeE = length(tempC(1,:,1));
                addon = [];
            else
                timeS = length(tempI(1,:,1))-length(tempC(1,:,1));
                timeE = length(tempI(1,:,1));
                addon = NaN(length(tempC(:,1,1)),timeS,length(tempC(1,1,:)));
            end                        
            temp = tempC;
            tempC = cat(2,addon,temp);
        end
    else
        timeS = 1;
        timeE = length(tempC(1,:,1));
    end
    out.MeanCorrect = NaN(length(tempC(1,1,:)),timeE-timeS+1);
    out.MeanIncorrect = NaN(length(tempC(1,1,:)),timeE-timeS+1);
    out.ConfidenceIntervalCorrect = NaN(length(tempC(1,1,:)),timeE-timeS+1);
    out.ConfidenceIntervalIncorrect = NaN(length(tempC(1,1,:)),timeE-timeS+1);
    %Individual cell level
%     f = figure(1);
%     f.Visible = 'off';
    for cellNum = 1 : length(tempC(1,1,:))
        if ~isempty(find(indpassed == cellNum,1))
%             figure(1)
            %Calculate mean Fluorescence and Confidence Intervals
            for time = timeS : timeE
                plotCtemp = tempC(:,time,cellNum);
                plotC(time) = nanmean(plotCtemp,1);
                confIntC(time) = (std(tempC(:,time,cellNum))/(sqrt(length(tempC(:,1,1)))))*1.96;
                if ~isempty(tempI)
                    plotItemp = tempI(:,time,cellNum);
                    plotI(time) = nanmean(plotItemp,1);
                    confIntI(time) = (std(tempI(:,time,cellNum))/(sqrt(length(tempI(:,1,1)))))*1.96;
                end
            end
            %Generate figures
            p1= shadedErrorBar((timeS:time),plotC(timeS:time),confIntC(timeS:time),'lineprops', '-g');
            p1.edge(1).HandleVisibility = 'off';
            p1.edge(2).HandleVisibility = 'off';
            hold on            
            if ~isempty(tempI)
                p2 = shadedErrorBar((timeS:time),plotI(timeS:time),confIntI(timeS:time),'lineprops', '-r');
                p2.edge(1).HandleVisibility = 'off';
                p2.edge(2).HandleVisibility = 'off';
                legend('Correct 95% Confidence Interval','Mean Correct Fluorescence','Incorrect 95% Confidence Interval','Mean Incorrect Fluorescence')
            else 
                legend('Correct 95% Confidence Interval','Mean Correct Fluorescence')
            end
            if cut == 1
                title(['Delay Period Confidence Intervals, Position: ', num2str(pos),' Cell: ', num2str(cellNum)])
            elseif cut == 2
                title(['Front Anchored Confidence Intervals, Position: ', num2str(pos),' Cell: ', num2str(cellNum)])
            elseif cut ==3
                title(['Back Anchored Confidence Intervals, Position: ', num2str(pos),' Cell: ', num2str(cellNum)])
            end
            xticks(timeS:fps:time);
            xticklabels(Frame2SecLabels(length(timeS:time),fps,ticmultiplier));            
            xlabel('Time (Seconds)')
            ylabel('df/f')            
            xlim([timeS time])
%             legend('Correct 95% Confidence Interval','Mean Correct Fluorescence','Incorrect 95% Confidence Interval','Mean Incorrect Fluorescence')
            saveas(gcf, [path,'/ConfidenceIntervals/',num2str(pos),'/',num2str(cellNum),'ConfidenceInterval.jpg'])
            out.MeanCorrect(cellNum,:) = plotC(timeS:timeE);
            if ~isempty(tempI)
                out.MeanIncorrect(cellNum,:) = plotI(timeS:timeE);
            else
                out.MeanIncorrect(cellNum,:) = NaN(1,timeE-timeS+1);
                confIntI = NaN(1,timeE-timeS+1);
            end
            out.ConfidenceIntervalCorrect(cellNum,:) = confIntC(timeS:timeE);
            out.ConfidenceIntervalIncorrect(cellNum,:) = confIntI(timeS:timeE);
            clf
        end
    end
    %Population level
    %Calculate mean Fluorescence and Confidence Intervals
    popCorrect = nanmean(out.MeanCorrect,1);
    popIncorrect = nanmean(out.MeanIncorrect,1);
    confIntCpop = [] ;
    confIntIpop = [];
    for time = 1 : length(out.MeanIncorrect(1,:))
        confIntCpop(time) = (std(out.MeanCorrect(indpassed,time))/(sqrt(length(out.MeanCorrect(indpassed,1)))))*1.96;
        confIntIpop(time) = (std(out.MeanIncorrect(indpassed,time))/(sqrt(length(out.MeanIncorrect(indpassed,1)))))*1.96;
    end
    %Generate figures
%     f = figure(1);
%     f.Visible = 'off';
    p1= shadedErrorBar((1:time),popCorrect,confIntCpop,'lineprops', '-g');
    p1.edge(1).HandleVisibility = 'off';
    p1.edge(2).HandleVisibility = 'off';
    hold on
    p2 = shadedErrorBar((1:time),popIncorrect,confIntIpop,'lineprops', '-r');
    p2.edge(1).HandleVisibility = 'off';
    p2.edge(2).HandleVisibility = 'off';    
    if cut == 1
        title(['Population Level Delay Period Confidence Intervals, Position: ', num2str(pos)])
    elseif cut == 2
        title(['Population Level Front Anchored Confidence Intervals, Position: ', num2str(pos)])
    elseif cut ==3
        title(['Population Level Back Anchored Confidence Intervals, Position: ', num2str(pos)])
    end
    xticks(timeS:fps:time);
    xticklabels(Frame2SecLabels(length(timeS:time),fps,ticmultiplier));            
    xlabel('Time (Seconds)')
    ylabel('df/f')    
    xlim([timeS time])
    legend('Correct 95% Confidence Interval','Population Mean Correct Fluorescence','Incorrect 95% Confidence Interval','Population Mean Incorrect Fluorescence')
    saveas(gcf, [path,'/ConfidenceIntervals/',num2str(pos),'/','POPULATIONConfidenceInterval.jpg'])
    out.populationMeanCorrect = popCorrect;
    out.populationMeanIncorrect = popIncorrect;
    out.populationConfidenceIntervalCorrect = confIntCpop;
    out.populationConfidenceIntervalIncorrect = confIntIpop;
    MatStruct = out;
    save([path,'/ConfidenceIntervals/',num2str(pos),'/MatStruc.mat'],'MatStruct')
    clf
end
end
