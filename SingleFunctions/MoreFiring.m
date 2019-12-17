function [ms] = MoreFiring(ms)
%%Produces firing activity of Calcium trace using Pnevmatikakis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Takes calcium trace from 'ms' file input and outputs its respective spike%
%activity in accordance with Pnevmatikakis-"Simultaneous Denoising,       %
%Deconvolution, and Demixing of Calcium Imaging Data". Need to add        %
%"contrained_foopsi.m" and "cvx_foopsi.m" to the path. These can be       %
%downloaded from "https://github.com/flatironinstitute/CaImAn-MATLAB".    %
%From here multiple methods should be able to be implimented.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author Emmanuel Wilson

trace = ms.trace;
firing_cnmfe = zeros(length(trace(:,1)), length(trace(1,:)));

for cellNum = 1 : length(trace(1,:))
    [~,s] = deconvolveCa(trace(:,cellNum));
    firing_cnmfe(:,cellNum) = s;
end
ms.firing_cnmfe = firing_cnmfe;
end