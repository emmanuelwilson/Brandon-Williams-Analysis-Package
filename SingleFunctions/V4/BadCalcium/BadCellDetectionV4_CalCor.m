function [CalCor] = BadCellDetectionV4_CalCor(ms)
CalCor = zeros(length(ms.FiltTraces(1,:)));
for i = 1 : length(ms.FiltTraces(1,:))
    for j = i + 1 : length(ms.FiltTraces(1,:))
        CalCor(i,j) = corr2(ms.FiltTraces(:,i),ms.FiltTraces(:,j));
    end
end
end