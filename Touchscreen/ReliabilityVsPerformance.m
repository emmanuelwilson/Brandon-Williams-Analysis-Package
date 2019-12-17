%% 6) Single cell Reliability as a function of task performance
%% Step 1) For each cell registered across days
folderpath = pwd;
load('aligmentMaps','Singlemap')
bmap = Singlemap;
bmap(bmap > 0) = 1;
elim = find(sum(bmap,2)==1);
ind = 1 : length(Singlemap(:,1));
ind(elim) = [];

% MFdelay = MeanFluoAcrossDays(folderpath,'Delay',1);
% MFfront = MeanFluoAcrossDays(folderpath,'Front',2);
MFback = MeanFluoAcrossDays(folderpath,'Back',3);
Rely = ReliabilityAcrossDays(folderpath);

mkdir([folderpath,'/Reliability/Delay'],'PerformanceVsReliability');
mkdir([folderpath,'/Reliability/FrontAnchored'],'PerformanceVsReliability');
mkdir([folderpath,'/Reliability/BackAnchored'],'PerformanceVsReliability');

dscat = NaN(length(Singlemap(1,:)),2);
dscattot = NaN(length(Singlemap(1,:)),2);
fscat = NaN(length(Singlemap(1,:)),2);

fscattot = NaN(length(Singlemap(1,:)),2);

bscat = NaN(length(Singlemap(1,:)),2);
bscattot = NaN(length(Singlemap(1,:)),2);
m = 5;
n = 4;
c = 1;

%Delay
for i = 1 : length(ind)    
    for j = 1 : length(Singlemap(1,:))        
        if ~isnan(Rely.dCorrectRely(ind(i),j))
            dscat(j,1) = Rely.dCorrectRely(ind(i),j);
            dscat(j,2) = MFdelay.BehavePerformance(1,j)/(MFdelay.BehavePerformance(1,j) + MFdelay.BehavePerformance(2,j));            
        end                
    end 
    
    if i == 1 
        dscattot = dscat;
    else
        dscattot = cat(1,dscattot,dscat);
    end
    
    figure(1)
    subplot(m,n,c)
    c = c + 1;
    scatter(dscat(:,1),dscat(:,2))
    ylim([0 1])
    xlim([0 1])
    xlabel('Reliability (splithalf Correlation)')
    ylabel('Performance')
    title(['Performance vs Reliability: Cell: ', num2str(ind(i))])
    
    if c > 20 
        c = 1; 
        saveas(gcf,[folderpath,'/Reliability/Delay/PerformanceVsReliability/Cell', num2str(ind(i)),'DelayPvsR.jpg']);
        pause(0.01)
        clf
    end
end

%Front Anchored
for i = 1 : length(ind)    
    for j = 1 : length(Singlemap(1,:))        
        if ~isnan(Rely.fCorrectRely(ind(i),j))
            fscat(j,1) = Rely.fCorrectRely(ind(i),j);
            fscat(j,2) = MFfront.BehavePerformance(1,j)/(MFfront.BehavePerformance(1,j) + MFfront.BehavePerformance(2,j));
        end                
    end
    
    if i == 1 
        fscattot = fscat;
    else
        fscattot = cat(1,fscattot,fscat);
    end
    
    figure(1)
    subplot(m,n,c)
    ylim([0 1])
    xlim([0 1])
    c = c + 1;
    scatter(fscat(:,1),fscat(:,2))
    xlabel('Reliability (splithalf Correlation)')
    ylabel('Performance')
    title(['Performance vs Reliability: Cell: ', num2str(ind(i))])
    
    if c > 20 
        c = 1; 
        saveas(gcf,[folderpath,'/Reliability/FrontAnchored/PerformanceVsReliability/Cell', num2str(ind(i)),'FrontPvsR.jpg']);
        pause(0.01)
        clf
    end
end

%Back Anchored
for i = 1 : length(ind)    
    for j = 1 : length(Singlemap(1,:))        
        if ~isnan(Rely.bCorrectRely(ind(i),j))
            bscat(j,1) = Rely.bCorrectRely(ind(i),j);
            bscat(j,2) = MFback.BehavePerformance(1,j)/(MFback.BehavePerformance(1,j) + MFback.BehavePerformance(2,j));
        end                
    end
    
    if i == 1 
        bscattot = bscat;
    else
        bscattot = cat(1,bscattot,bscat);
    end
    
    figure(1)
    subplot(m,n,c)
    c = c + 1;
    scatter(bscat(:,1),bscat(:,2))
    ylim([0 1])
    xlim([0 1])
    xlabel('Reliability (splithalf Correlation)')
    ylabel('Performance')
    title(['Performance vs Reliability: Cell: ', num2str(ind(i))])
    
    if c > 20 
        c = 1; 
        saveas(gcf,[folderpath,'/Reliability/BackAnchored/PerformanceVsReliability/Cell', num2str(ind(i)),'BackPvsR.jpg']);
        pause(0.01)
        clf
    end
end

