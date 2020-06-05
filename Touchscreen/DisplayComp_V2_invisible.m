%% Will Display and save Figures of cells with all trial types side by side
% INPUT: 
%   -rasts: structure of Cell arrays containing calcium. Each column in
%   cell array represents a choice made, each matrix in each cell NxMxC 
%   where N represents the trial, M for time in frames and C for the cell 
%   index
%   -minmax: minimum and maximum value for each cell calcium transient
%   -name: Folder name where results will be saved
%   -indpassed: index of cells which passed criteria
% OUTPUT:
%   -Saves figures in respective folders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Émmanuel Wilson

function [] = DisplayComp_V2_invisible(rasts,minmax, name,indpassed)

fps = 30;
ticmultiplier = 2;
maxval = 0;
m1 = 0;
m2 = 0;
m3 = 0;
m4 = 0;
mli = [];
m2i = [];
m3i = [];
m4i = [];
fnames = fieldnames(rasts);

for i = 1: length(rasts.CorrectTrial)    
    if any(strcmp(fnames, 'CorrectTrial'))
        m1 = m1 + ~isempty(rasts.CorrectTrial{i});
        m1i(i) = ~isempty(rasts.CorrectTrial{i});
    end
    if any(strcmp(fnames, 'IncorrectTrial'))
        m2 = m2 + ~isempty(rasts.IncorrectTrial{i});
        m2i(i) = ~isempty(rasts.IncorrectTrial{i});
    end
    if any(strcmp(fnames, 'CorrectCorrectionTrial'))
        m3 = m3 + ~isempty(rasts.CorrectCorrectionTrial{i});
        m3i(i) = ~isempty(rasts.CorrectCorrectionTrial{i});
    end
    if any(strcmp(fnames, 'IncorrectCorrectionTrial'))
        m4 = m4 + ~isempty(rasts.IncorrectCorrectionTrial{i});
        m4i(i) = ~isempty(rasts.IncorrectCorrectionTrial{i});
    end
end

maxval = max(cat(1,m1,m2,m3,m4));
if m1 == maxval
    poke = find(m1i);
elseif m2 == maxval
    poke = find(m2i);
elseif m3 == maxval
    poke = find(m3i);
elseif m4 == maxval
    poke = find(m4i);
end

cellnum = length(rasts.CorrectTrial{poke(1)}(1,1,:));
minlim = minmax(:,1);
maxlim = minmax(:,2);
n = maxval;
m = 4;
%Individual cell save
mkdir(name,'ComparisonView')
f = figure(1);
f.Visible = 'off';
for i = 1: cellnum
    if ~isempty(find(indpassed == i,1))
    c = 1;
    for j = 1: maxval
        pos = poke(j);
        f = figure(1);
%         f.Visible = 'off';
%         colormap parula
        subplot(n,m,c)
        c = c + 1;
        if ~isempty(rasts.CorrectTrial{pos})
            imagesc(rasts.CorrectTrial{pos}(:,:,i))
            xticks(1:fps:length(rasts.CorrectTrial{pos}(:,:,i)));
            xticklabels(Frame2SecLabels(length(rasts.CorrectTrial{pos}(:,:,i)),30,ticmultiplier));
            yticks(1:2:length(rasts.CorrectTrial{pos}(:,:,i)));
        end
        colorbar
        caxis([minlim(i) maxlim(i)])
        title(['Correct Trials ', num2str(pos)])
        
        xlabel('Time(Seconds)')
        ylabel('Trial')
        
        subplot(n,m,c)
        c = c + 1;
        if ~isempty(rasts.IncorrectTrial{pos})
            imagesc(rasts.IncorrectTrial{pos}(:,:,i))
            xticks(1:fps:length(rasts.IncorrectTrial{pos}(:,:,i)));
            xticklabels(Frame2SecLabels(length(rasts.IncorrectTrial{pos}(:,:,i)),30,ticmultiplier));
            yticks(1:length(rasts.IncorrectTrial{pos}(:,:,i)));
        end
        colorbar
        caxis([minlim(i) maxlim(i)])
        title(['Incorrect Trials ', num2str(pos)])        
        xlabel('Time(Seconds)')
        ylabel('Trial')
        
        subplot(n,m,c)
        c = c + 1;
        if ~isempty(rasts.CorrectCorrectionTrial{pos})
            imagesc(rasts.CorrectCorrectionTrial{pos}(:,:,i))
            xticks(1:fps:length(rasts.CorrectCorrectionTrial{pos}(:,:,i)));
            xticklabels(Frame2SecLabels(length(rasts.CorrectCorrectionTrial{pos}(:,:,i)),30,ticmultiplier));
            yticks(1:length(rasts.CorrectCorrectionTrial{pos}(:,:,i)));
        end
        colorbar
        caxis([minlim(i) maxlim(i)])
        title(['Correct Correction Trials ', num2str(pos)])        
        xlabel('Time(Seconds)')
        ylabel('Trial')
        
        subplot(n,m,c)
        c = c + 1;
        if ~isempty(m4i) && ~isempty(rasts.IncorrectCorrectionTrial{pos})
            imagesc(rasts.IncorrectCorrectionTrial{pos}(:,:,i))
            xticks(1:fps:length(rasts.IncorrectCorrectionTrial{pos}(:,:,i)));
            xticklabels(Frame2SecLabels(length(rasts.IncorrectCorrectionTrial{pos}(:,:,i)),30,ticmultiplier));
            yticks(1:length(rasts.IncorrectCorrectionTrial{pos}(:,:,i)));
        end
        colorbar
        caxis([minlim(i) maxlim(i)])
        title(['Incorrect Correction Trials ', num2str(pos)])        
        xlabel('Time(Seconds)')
        ylabel('Trial')
    end
    saveas(gcf,[name ,'/ComparisonView/Cell',num2str(i),'.jpeg' ]);
    end
end

%population level figures
c = 1;
f = figure(1);
f.Visible = 'off';
for i = 1 : maxval
    pos = poke(i);
    if ~isempty(rasts.CorrectTrial{pos})
        cor = mean(rasts.CorrectTrial{pos}(:,:,indpassed),1);
        cor = permute(cor,[3 2 1]);
        cor = sortpeaks(cor);
    else
        cor = [];
    end    
    if ~isempty(rasts.IncorrectTrial{pos})
        icor = mean(rasts.IncorrectTrial{pos}(:,:,indpassed),1);
        icor = permute(icor,[3 2 1]);
        icor = sortpeaks(icor);
    else
        icor = [];
    end
    if ~isempty(rasts.CorrectCorrectionTrial{pos})
        ccor = mean(rasts.CorrectCorrectionTrial{pos}(:,:,indpassed),1);
        ccor = permute(ccor,[3 2 1]);
        ccor = sortpeaks(ccor);
    else 
        ccor = [];
    end
    if ~isempty(m4i) && ~isempty(rasts.IncorrectCorrectionTrial{pos})
        iccor = mean(rasts.IncorrectCorrectionTrial{pos}(:,:,indpassed),1);
        iccor = permute(iccor,[3 2 1]);
        iccor = sortpeaks(iccor);
    else
        iccor = [];
    end
    
%     figure(1)
    colormap parula
    subplot(n,m,c)
    c = c + 1;
    imagesc(cor(:,:))
    colorbar
    title(['Correct Trials ', num2str(pos)])
    if ~isempty(cor)
        xticks(1:fps:length(cor(1,:)));
        xticklabels(Frame2SecLabels(length(cor(1,:)),30,ticmultiplier));
        yticks(1:round(length(cor(:,1))/10):length(cor(:,1)));
    end
    xlabel('Time(Seconds)')
    ylabel('Cell ID')
    
    subplot(n,m,c)
    c = c + 1;
    imagesc(icor(:,:))
    colorbar
    title(['Incorrect Trials ', num2str(pos)])
    if ~isempty(icor)
        xticks(1:fps:length(icor(1,:)));
        xticklabels(Frame2SecLabels(length(icor(1,:)),30,ticmultiplier));
        yticks(1:round(length(icor(:,1))/10):length(icor(:,1)));
    end
    xlabel('Time(Seconds)')
    ylabel('Cell ID')
    
    subplot(n,m,c)
    c = c + 1;
    imagesc(ccor(:,:))
    colorbar
    title(['Correct Correction Trials ', num2str(pos)])
    if ~isempty(ccor)
        xticks(1:fps:length(ccor(1,:)));
        xticklabels(Frame2SecLabels(length(ccor(1,:)),30,ticmultiplier));
        yticks(1:round(length(ccor(:,1))/10):length(ccor(:,1)));
    end
    xlabel('Time(Seconds)')
    ylabel('Cell ID')
    
    subplot(n,m,c)
    c = c + 1;
    imagesc(iccor(:,:))
    colorbar
    title(['Incorrect Correction Trials ', num2str(pos)])
    if ~isempty(icor)
        xticks(1:fps:length(icor(1,:)));
        xticklabels(Frame2SecLabels(length(icor(1,:)),30,ticmultiplier));
        yticks(1:round(length(icor(:,1))/10):length(icor(:,1)));
    end
    xlabel('Time(Seconds)')
    ylabel('Cell ID')
    
end

saveas(gcf,[name ,'/ComparisonView/PopulationLevel.jpeg' ]);

end

function rast = sortpeaks(rast)
[~,maxind] = max(rast,[],2);
[~, rastsort] = sort(maxind);
rast = rast(rastsort,:);
end