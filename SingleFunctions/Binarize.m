%%This function used the data saved in ms for a cell firing and converts the transients to binary. 
%Author: Tori-Lynn Temple 
    
function  spatialCoding = Binarize(ms) 
 
% dt = median(diff(ms.time))/1000; % Conversion from ms to s
% Fs = 1/dt;
z_threshold = 2;                                                            %Zscore Threshold, typically ~2

[bFilt,aFilt] = butter(2,  2/(30/2), 'low');                                %LowPassFilter parameter creation

for trace_i = 1: length(ms.FiltTraces(1,:))                                 %Cycle through all cells
raw_trace = ms.FiltTraces(:,trace_i);                                       %select trace
filt_trace = zscore(filtfilt(bFilt,aFilt,raw_trace));                       %Filter trace and Zscore result
d1_trace = diff(filt_trace);                                                %Derivative of trace
d1_trace(end+1) = 0;
d2_trace = diff(d1_trace);                                                  %Derivative of derived trace
d2_trace(end+1) = 0;

binary_trace = filt_trace*0;                                                %initialize Binarized Trace
binary_trace(filt_trace>z_threshold & d1_trace>0) = 1;                      %Set all positive slope transients above the z threshold to 1

spatialCoding.binarizedTraces(:, trace_i) = binary_trace;                   %output

end 