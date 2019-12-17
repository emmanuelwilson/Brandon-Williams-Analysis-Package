%% Smooth deconvolved trace using Gaussian filter where future values don't effect past values.
%INPUT:
%   -deconvolvedTraces: n by m matrix of the deconvolved calcium traces
%OUTPUT:
%   -smoothedTrace: n by m matrix of the smoothed deconvolvedTraces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Author: Émmanuel Wilson

function smoothedTraces = HalfGaussFilt(deconvolvedTraces)
x = [-1.5:0.1:1.5];
fullgauss = normpdf(x,0,0.5);
halfgauss = fullgauss(round(length(fullgauss)/2):end);

for i =1 : length(deconvolvedTraces(1,:))
    smoothedTraces(:,i) = conv(deconvolvedTraces(:,i),halfgauss);
end

end