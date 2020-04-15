function [v gapSize] = interpNaNs(v)
    gapSize = zeros(length(v(:,1)),1);
    while any(isnan(v(:)))
        start = find(any(isnan(v(:,1)),2),1,'first');
        
        % Solve End Edge Case
        if start==length(v(:,1))
            stop = 1;
        else
            stop = find(~any(isnan(v(start+1:end,:)),2),1,'first');
        end
        
        % Solve End Edge Case
        if isempty(stop) || length(v(:,1)) == start+stop-1
            for i = start:length(v(:,1))
                v(i,:) = v(start-1,:);
            end
            gapSize(start:length(v(:,1))) = length(v(:,1)) - start;
            continue
        end
        
        % Solve Beginning Edge Case
        if start == 1
            for i = start:start+stop-1
                v(i,:) = v(start+stop,:);
            end
            gapSize(start:start+stop-1) = start+stop-1;
            continue
        end
        
        % Linterp missing values
        for i = 1:length(v(1,:))
            a = linterp([1 stop+2],[v(start-1,i) v(start+stop,i)],[1:stop+2]);
            v(start:start+stop-1,i) = a(2:end-1);
            gapSize(start:start+stop-1) = stop;
        end
    end
end