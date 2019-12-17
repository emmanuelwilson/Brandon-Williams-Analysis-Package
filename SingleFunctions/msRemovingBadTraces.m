%Author: Tori-Lynn Temple, Janurary 31st, 2018

function [ms] = msRemovingBadTraces(ms)
%this function runs after found traces using ____. From this we use a
%simple algorithm to identify noise traces in individual cells and remove
%them from the analysis.  

current_trace = ms.trace(:,trace_i); % uploads the current trace of a cell  
                                     % into array current trace 

%this goes through a trace to find all the peaks that are atleast 1 __ high
%below the largest peak. 
[peaks, locationOfPeaks] =  findpeaks(current_trace, Fs, 'Threshold' 1); 
 
if isempty(peaks) 
    % this is where the traces are deleted if resulting in uneven trace data
    % Meaning the traces are insufficient for cell characteristics. The cell
    % for which the trace belongs too will then be deleted, it is most likely 
    % not a cell and can be neglected from further analysis. 
    rmfield(ms.transients{trace_i} = [];)
   
end 

end 