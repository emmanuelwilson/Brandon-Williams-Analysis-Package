%% Correlates correct trial with subsampled correct trial (control metric)
% INPUT:
%   -rasts: Raster plots of all decisions made
% OUTPUT:
%   -out: N x M matrix, where N is the cell# and M the mean correlation
%   value of the subsampled calcium for each choice. ie, if only two
%   choices were made (1 or 5) then the first column is the mean 
%   correlation for choice1 and the second column the mean correlation for 
%   choice5. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Émmanuel Wilson

function out = tMetric2(rasts)

maxval = 0;
m1 = 0;
m2 = 0;

for i = 1: length(rasts.CorrectTrial)
    m1 = m1 + ~isempty(rasts.CorrectTrial{i});
    m1i(i) = ~isempty(rasts.CorrectTrial{i});
    m2 = m2 + ~isempty(rasts.IncorrectTrial{i});
    m2i(i) = ~isempty(rasts.IncorrectTrial{i});
end

maxval = max(cat(1,m1,m2)); %Max number of choices made

%Find which trial had the most choices
if m1 == maxval
    poke = find(m1i);
elseif m2 == maxval
    poke = find(m2i);
end

out = zeros(length(rasts.CorrectTrial{poke(1)}(1,1,:)),maxval,100);

%Correlation between every population raster plot
for i = 1: maxval
    pos1 = poke(i);
    tempC = ExtractTrial(rasts,1,pos1);
    tempI = ExtractTrial(rasts,2,pos1);
    if ~isempty(tempI) && ~isempty(tempC)
        nI = length(tempI(:,1,1));
        nC = length(tempC(:,1,1));
        if nI< nC
            for j = 1 : 100
                randTrials = randperm(nC,nI);
                comp1 = tempC;
                comp1(randTrials,:,:) = [];
                comp2 = tempC(randTrials,:,:);
                comp1 = nanmean(comp1,1);
                comp1 = permute(comp1,[3 2 1]);
                comp2 = nanmean(comp2,1);
                comp2 = permute(comp2, [3 2 1]);
                for cellNum = 1: length(comp1(:,1))
                    out(cellNum,i,j) = corr2(comp1(cellNum,:),comp2(cellNum,:));
                    if (isnan(out(cellNum,i,j)))
                        out(cellNum,i,j) = 0;
                    end
                end
            end
        end
    end
    out = mean(out,3);
end
end