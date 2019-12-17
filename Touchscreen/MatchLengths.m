%Will match the length of all trials with the shortest trial length,
%dropping any extra components
function [out] = MatchLengths(trial,cl,cut,ind,dind,dframe)
if ~isempty(trial)
    if ~(length(trial(1,1,:)) == 0)
        out = zeros(size(trial));
        if ~isempty(trial) && length(trial(1,1,:))>1
            temp = permute(trial,[3 2 1]);
            if cut == 2 %Front Anchored: eliminate back end excess
                if ~(min(temp) == length(temp))
                    out = permute(temp(:,1:min(cl),:),[3 2 1]);
                end
            elseif cut == 3 && ~isempty(ind) && ~isempty(dind) && ~isempty(dframe) %Back Anchored: eliminate front end excess
                sT = cl - ((dind(:) - ind(:))+dframe);                      %Shortest Trial for the BA is found after subtracting the delay period
                if ~(min(temp) == length(temp))
                    out = permute(temp(:,end-min(sT)+1:end,:),[3 2 1]);
                end
            end            
        end
    else
        out = NaN;
    end
else
    out = NaN;
end    
end