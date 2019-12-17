
function cormat = cormatrix(folderpath)

folder = dir(folderpath);
n = length(folder(:,1))-2;
tn = roots([1 -1 -n*2]);
trialNum = tn(1);
d = ones(1,int32(trialNum));
cormat = diag(d);

for i = 1 : trialNum-1
    for j = i+1 : trialNum
%         loc = find(contains(folder.name,[num2str(i),'_',num2str(j)]));
        corval = load([folderpath,'\',num2str(i),'_',num2str(j),'\correlation.mat']);
        cormat(i,j) = corval.maximal_cross_correlation;
        cormat(j,i) = corval.maximal_cross_correlation;
    end
end

end