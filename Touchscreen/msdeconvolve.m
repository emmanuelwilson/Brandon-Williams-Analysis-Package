%% Adds Second Order auto-Regression signal and denoised Ca trace to the ms structure
%INPUT:
%   -ms: miniscope data structure, must contain ms.numNeurons and
%   ms.FiltTraces. 
%OUTPUT:
%   -ms: same strucutre with added feilds, denoisedCa and deconvolvedSig.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Émmanuel Wilson

function [ms] = msdeconvolve(ms)

[ms.deconvolution.denoisedCa(:,1),ms.deconvolution.deconvolvedSig(:,1),ms.deconvolution.options] = deconvolveCa(ms.FiltTraces(:,1), 'ar1');

for i = 2 : ms.numNeurons
    [ms.deconvolution.denoisedCa(:,i),ms.deconvolution.deconvolvedSig(:,i)] = deconvolveCa( ms.FiltTraces(:,i), 'ar1');
end

end