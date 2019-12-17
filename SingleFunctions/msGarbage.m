%Eliminates specified cells from footprints and miniscope structure
function [ms, SFP] = msGarbage(ms,SFP,elim)

ms.FiltTraces(:,elim) = [];
ms.RawTraces(:,elim) = [];
ms.SFPs(:,:,elim) = [];
ms.numNeurons = ms.numNeurons - length(elim);
SFP(elim,:,:) = [];

end