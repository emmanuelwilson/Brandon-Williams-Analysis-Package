%Looks if folder exists, otherwise makes folder
function [] = checkP(path)

if ~exist(path)
    mkdir(path);    
end

end