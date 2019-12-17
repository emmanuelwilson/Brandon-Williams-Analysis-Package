function [out] = SequenceViewer(Calcium,Cellpassed,minmax)
sequence = [];
neuron = [];
for i = 1 : length(Cellpassed)    
    for j = 1 : length(Calcium(1,1,:))
        if i == 1 && j == 1 && ~isempty(find(isnan(Calcium(Cellpassed(i),:,1)),1))
            neuron = Calcium(Cellpassed(i),1:find(isnan(Calcium(Cellpassed(i),:,j)),1),j);
        elseif i == 1 && j == 1
            neuron = Calcium(Cellpassed(i),:,j);
        elseif ~isempty(find(isnan(Calcium(Cellpassed(i),:,j)),1))
            seq = Calcium(Cellpassed(i),1:find(isnan(Calcium(Cellpassed(i),:,j)),1),j);
            neuron = cat(2,neuron,seq);
        else
            seq = Calcium(Cellpassed(i),:,j);
            neuron = cat(2,neuron,seq);
        end
    end        
    if i == 1
        sequence = neuron;
    else
        sequence = cat(1,sequence,neuron);
    end
    neuron = [];
end
sequence = (sequence - minmax(Cellpassed,1)) ./ (minmax(Cellpassed,2) - minmax(Cellpassed,1));
out = sequence;
end