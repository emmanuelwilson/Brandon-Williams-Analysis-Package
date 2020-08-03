%% Adds Second Order auto-Regression signal and denoised Ca trace to the ms structure
%INPUT:
%   -ms: miniscope data structure, must contain ms.numNeurons and
%   ms.FiltTraces. 
%OUTPUT:
%   -ms: same strucutre with added feilds, denoisedCa and deconvolvedSig.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Émmanuel Wilson

function [ms] = msdeconvolve(ms)

for i = 1 : ms.numNeurons
    [ms.deconvolution.denoisedCa(:,i),ms.deconvolution.deconvolvedSig(:,i)] = deconvolveCa( ms.FiltTraces(:,i), 'ar2');
end

end