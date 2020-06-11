%% Pie chart generation for TUNL task, generate numbers for all combinations of cells which pass Delay and Reward criteria for Seperation 3
%INPUT: Choice seperated Shuffled criteria
%Output: 
%   R1 : R = Reward, 1 = choice 1
%   D1 : D = Delay, 1 = choice 1

function [R1,R5,D1,D5,R1D1,R5D1,R5D5,R1D5,R1R5,D1D5,R1R5D1,R1R5D5,R1D1D5,R5D1D5,R1R5D1D5,totalCells] = S3_PieChartVariables(shuffledCritChoice)

R1 = shuffledCritChoice.BcorrectCells{1};
R5 = shuffledCritChoice.BcorrectCells{5};
D1 = shuffledCritChoice.DcorrectCells{1};
D5 = shuffledCritChoice.DcorrectCells{5};
totalCells = legnth(shuffledCritChoice.Bpercentdic(:,1));

%find all common values
R1D1 = intersect(R1,D1);
R1D5 = intersect(R1,D5);
R5D1 = intersect(R5,D1);
R5D5 = intersect(R5,D5);
R1R5 = intersect(R1,R5);
D1D5 = intersect(D1,D5);
R1R5D1 = intersect(R1R5,D1);
R1R5D5 = intersect(R1R5,D5);
R1D1D5 = intersect(R1,D1D5);
R5D1D5 = intersect(R5,D1D5);
R1R5D1D5 = intersect(R1R5,D1D5);

%eliminate cells which appear in multiple variables
if ~isemtpy(R1R5D1D5)    
    [~,R1ind,~] = intersect(R1,R1R5D1D5);
    R1(R1ind) = [];
    [~,R5ind,~] = intersect(R5,R1R5D1D5);
    R5(R5ind) = [];
    [~,D1ind,~] = intersect(D1,R1R5D1D5);
    D1(D1ind) = [];
    [~,D5ind,~] = intersect(D5,R1R5D1D5);
    D5(D5ind) = [];
    [~,R1R5ind,~] = intersect(R1R5,R1R5D1D5);
    R1R5(R1R5ind) = [];
    [~,D1D5ind,~] = intersect(D1D5,R1R5D1D5);
    D1D5(D1D5ind) = [];    
    [~,R1D1ind,~] = intersect(R1D1,R1R5D1D5);
    R1D1(R1D1ind) = [];
    [~,R5D1ind,~] = intersect(R5D1,R1R5D1D5);
    R5D1(R5D1ind) = [];
    [~,R1D5ind,~] = intersect(R1D5,R1R5D1D5);
    R1D5(R1D5ind) = [];
    [~,R5D5ind,~] = intersect(R5D5,R1R5D1D5);
    R5D5(R5D5ind) = [];
    [~,R1R5D1ind,~] = intersect(R1R5D1,R1R5D1D5);
    R1R5D1(R1R5D1ind) = [];
    [~,R1R5D5ind,~] = intersect(R1R5D5,R1R5D1D5);
    R1R5D5(R1R5D5ind) = [];
    [~,R1D1D5ind,~] = intersect(R1D1D5,R1R5D1D5);
    R1D1D5(R1D1D5ind) = [];
    [~,R5D1D5ind,~] = intersect(R5D1D5,R1R5D1D5);
    R5D1D5(R5D1D5ind) = [];
    
end
if ~isemtpy(R1R5D1)    
    [~,R1ind,~] = intersect(R1,R1R5D1);
    R1(R1ind) = [];
    [~,R5ind,~] = intersect(R5,R1R5D1);
    R5(R5ind) = [];
    [~,D1ind,~] = intersect(D1,R1R5D1);
    D1(D1ind) = [];
    [~,R1R5ind,~] = intersect(R1R5,R1R5D1);
    R1R5(R1R5ind) = [];
end
if ~isemtpy(R1R5D5)    
    [~,R1ind,~] = intersect(R1,R1R5D5);
    R1(R1ind) = [];
    [~,R5ind,~] = intersect(R5,R1R5D5);
    R5(R5ind) = [];
    [~,D5ind,~] = intersect(D5,R1R5D5);
    D1(D5ind) = [];
    [~,R1R5ind,~] = intersect(R1R5,R1R5D5);
    R1R5(R1R5ind) = [];
end
if ~isemtpy(R1D1D5)    
    [~,R1ind,~] = intersect(R1,R1D1D5);
    R1(R1ind) = [];
    [~,D1ind,~] = intersect(D1,R1D1D5);
    D1(D1ind) = [];
    [~,D5ind,~] = intersect(D5,R1D1D5);
    D5(D5ind) = [];
    [~,D1D5ind,~] = intersect(D1D5,R1D1D5);
    D1D5(D1D5ind) = [];
end
if ~isemtpy(R5D1D5)    
    [~,R5ind,~] = intersect(R5,R5D1D5);
    R5(R5ind) = [];
    [~,D1ind,~] = intersect(D1,R5D1D5);
    D1(D1ind) = [];
    [~,D5ind,~] = intersect(D5,R5D1D5);
    D5(D5ind) = [];
    [~,D1D5ind,~] = intersect(D1D5,R5D1D5);
    D1D5(D1D5ind) = [];
end
if ~isemtpy(R1R5)    
    [~,R1ind,~] = intersect(R1,R1R5);
    R1(R1ind) = [];
    [~,R5ind,~] = intersect(R5,R1R5);
    R5(R5ind) = [];        
end
if ~isemtpy(D1D5)
    [~,D1ind,~] = intersect(D1,D1D5);
    D1(D1ind) = [];
    [~,D5ind,~] = intersect(D5,D1D5);
    D5(D5ind) = [];
end

if ~isemtpy(R1D1)    
    [~,R1ind,~] = intersect(R1,R1D1);
    R1(R1ind) = [];
    [~,D1ind,~] = intersect(D1,R1D1);
    D1(D1ind) = [];    
end

if ~isemtpy(R5D1)    
    [~,R5ind,~] = intersect(R5,R5D1);
    R5(R5ind) = [];
    [~,D1ind,~] = intersect(D1,R5D1);
    D1(D1ind) = [];    
end

if ~isemtpy(R1D5)
    [~,R1ind,~] = intersect(R1,R1D5);
    R1(R1ind) = [];
    [~,D5ind,~] = intersect(D5,R1D5);
    D5(D5ind) = [];    
end

if ~isemtpy(R5D5)
    [~,R5ind,~] = intersect(R5,R5D5);
    R5(R5ind) = [];
    [~,D5ind,~] = intersect(D5,R5D5);
    D5(D5ind) = [];    
end

end