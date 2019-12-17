%% Metric of correlation between correct decisions and all other decisions
%INPUT:
%   -rasts: Raster plots of all decisions made
%   -cut:
%       *1: Delay anchor
%       *2: FrontAnchor
%       *3: BackAnchor
% OUTPUT:
%   -out: N x N Cell array where N is the number of decision possibilities.
%   Each cell contains a structure which have the following variables: 
%   Correct, Incorrect, CorrectCorrection and IncorrectCorrection. In each 
%   subvariable contains a correlation value of correct Row# vs Column# choice
%   for each cell.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Émmanuel Wilson

function out = tMetric1(rasts,cut)
maxval = 0;
m1 = 0;
m2 = 0;
m3 = 0;
m4 = 0;
mli = [];
m2i = [];
m3i = [];
m4i = [];
ttype1 = 2:4;
ttype2 = 1:4;
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

if isempty(m2i)
    ttype1 = [];
    ttype2 = 1;
elseif isempty(m3i)
    ttype1 = ttype1([1,3]);
    ttype2 = ttype2([1,2,4]);
elseif isempty(m4i)
    ttype1 = ttype1(1:2);
    ttype2 = ttype2(1:3);
end

maxval = max(cat(1,m1,m2,m3,m4)); %Max number of choices made
out = cell(length(rasts.CorrectTrial),length(rasts.CorrectTrial));
%Find which trial had the most choices
if m1 == maxval
    poke = find(m1i);
elseif m2 == maxval
    poke = find(m2i);
elseif m3 == maxval
    poke = find(m3i);
elseif m4 == maxval
    poke = find(m4i);
end

%Correlation between every population raster plot
c1 = false;
c2 = false;
corrvals = cell(5,1);
for i = 1: maxval
    pos1 = poke(i);    
    for j = 1: maxval        
        pos2 = poke(j);
        if pos1 == pos2
            for t = ttype1
                comp1 = ExtractTrial(rasts,1,pos1);
                comp1 = permute(nanmean(comp1,1),[2 3 1]);
                comp2 = ExtractTrial(rasts,t,pos2);
                comp2 = permute(nanmean(comp2,1),[2 3 1]);
                if length(comp1(:,1)) > length(comp2(:,1))
                    c1 = true;
                elseif length(comp1(:,1)) < length(comp2(:,1))
                    c2 = true;
                end
                if cut == 2
                    if c1
                        comp1 = comp1(1:end-(length(comp1(:,1)) - length(comp2(:,1))),:);
                    elseif c2
                        comp2 = comp2(1:end-(length(comp2(:,1)) - length(comp1(:,1))),:);
                    end                    
                elseif cut == 3
                    if c1
                        comp1 = comp1((length(comp1(:,1)) - length(comp2(:,1))+1):end,:);
                    elseif c2
                        comp2 = comp2((length(comp2(:,1)) - length(comp1(:,1))+1):end,:);
                    end
                end
                c1 = false;
                c2 = false;
                if ~isempty(comp1) && ~isempty(comp2)
                    for c = 1 :length(comp1(1,:))
                        comp(c) = corr2(comp1(:,c),comp2(:,c));
                    end
                    corrvals{t} = comp;
                end                
            end            
        else
            for t = ttype2
                comp2 = ExtractTrial(rasts,t,pos2);
                comp2 = permute(nanmean(comp2,1),[2 3 1]);
                if length(comp1(:,1)) > length(comp2(:,1))
                    c1 = true;
                elseif length(comp1(:,1)) < length(comp2(:,1))
                    c2 = true;
                end
                if cut == 2
                    if c1
                        comp1 = comp1(1:end-(length(comp1(:,1)) - length(comp2(:,1))),:);
                    elseif c2
                        comp2 = comp2(1:end-(length(comp2(:,1)) - length(comp1(:,1))),:);
                    end                    
                elseif cut == 3
                    if c1
                        comp1 = comp1((length(comp1(:,1)) - length(comp2(:,1))+1):end,:);
                    elseif c2
                        comp2 = comp2((length(comp2(:,1)) - length(comp1(:,1))+1):end,:);
                    end
                end
                c1 = false;
                c2 = false;
                if ~isempty(comp1) && ~isempty(comp2)
                    for c = 1 :length(comp1(1,:))
                        comp(c) = corr2(comp1(:,c),comp2(:,c));
                    end
                    corrvals{t} = comp;
                end
            end
        end                    
        Correlations.Correct = corrvals{1};
        Correlations.Incorrect = corrvals{2};
        Correlations.CorrectCorrection = corrvals{3};
        Correlations.IncorrectCorrection = corrvals{4};
        
        out{pos1,pos2} = Correlations;
    end
end

end


% function [out,name] = ExtractTrial(rasts,num,choice)
% if num == 1
%     out = rasts.CorrectTrial{choice};    
% elseif num == 2
%     out = rasts.IncorrectTrial{choice};
% elseif num ==3
%     out = rasts.CorrectCorrectionTrial{choice};
% elseif num == 4
%     out = rasts.IncorrectCorrectionTrial{choice};
% end
% end