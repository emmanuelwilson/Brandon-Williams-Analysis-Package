%% Adds Second Order auto-Regression signal and denoised Ca trace to the ms structure
%INPUT:
%   -ms: miniscope data structure, must contain ms.numNeurons and
%   ms.FiltTraces. 
%OUTPUT:
%   -ms: same strucutre with added feilds, denoisedCa and deconvolvedSig.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Émmanuel Wilson

function [] = msdeconvolve_CC(p)

paths = genpath(p);
if ispc
    folders = strsplit(paths,';')';
else
    folders = strsplit(paths,':')';
end
for i = 1 : length(folders)
    if ~isempty(folders{i})
        d = dir(folders{i});
        fnames = {d.name}
        if ~isempty(find(strcmp(fnames,'ms.mat'),6)) %&& ~isempty(find(contains(folders{i},'Miniscope'),1))
            cd(folders{i})                                  %Change current folder
            try
                load([pwd '/ms.mat']);
                msnames = fieldnames(ms);
                if sum(contains(msnames,'deconvolvedSig')) == 0
                    x = [1:length(ms.FiltTraces(:,1))];
                    xq = [1:5:length(ms.FiltTraces(:,1))];
                    temp1 = interp1(x,ms.FiltTraces,xq);
                    ms.FiltTraces = interp1(xq,temp1,x);
                    ms.FiltTraces(find(isnan(ms.FiltTraces))) = 0;
                    ms = msdeconvolve(ms);
                    save([pwd '/ms.mat'],'ms','-v7.3');
                end
            catch
                display([folders{i},' Failed to analize'])
            end
        end
    end
end
end