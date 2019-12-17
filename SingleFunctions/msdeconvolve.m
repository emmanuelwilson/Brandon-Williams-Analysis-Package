%% Adds Second Order auto-Regression signal and denoised Ca trace to the ms structure
%INPUT:
%   -ms: miniscope data structure, must contain ms.numNeurons and
%   ms.FiltTraces. 
%OUTPUT:
%   -ms: same strucutre with added feilds, denoisedCa and deconvolvedSig.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: �mmanuel Wilson

function [ms] = msdeconvolve(ms)
for i = 1 : ms.numNeurons
    [ms.denoisedCa(:,i),ms.deconvolvedSig(:,i)] = deconvolveCa(ms.FiltTraces(:,i),'ar2',[1.416;-0.421776879189149]);
end
end