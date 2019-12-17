%% Mean direction of EBC ratemap
function [angle, angle1, angle2] = MRLangle(rm,rm1,rm2,deconvolved)
if ~isempty(rm)
    for i = 1 : length(rm(1,1,:))
        metric = mean(rm(:,:,i),1)';
        if ~deconvolved
            metric = metric - min(metric);
        end
        xs = metric(1:end-1).*cos(linspace(-180,179,length(rm(1,:,1))-1)); % average
        ys = metric(1:end-1).*sin(linspace(-180,179,length(rm(1,:,1))-1));
        
        angle = atan2(mean(ys),mean(xs)); % mean direction
    end
end

if ~isempty(rm1)
    for i = 1 : length(rm1(1,1,:))
        metric1 = mean(rm1(:,:,i),1)';
        if ~deconvolved
            metric1 = metric1 - min(metric1);
        end
        
        xs = metric1(1:end-1).*cos(linspace(-180,179,length(rm1(1,:,1))-1)); % average
        ys = metric1(1:end-1).*sin(linspace(-180,179,length(rm1(1,:,1))-1));
        
        angle1 = atan2(mean(ys),mean(xs)); % mean direction
    end
end

if ~isempty(rm2)
    for i = 1 : length(rm2(1,1,:))
        metric2 = mean(rm2(:,:,i),1)';
        if ~deconvolved
            metric2 = metric2 - min(metric2);
        end
        
        xs = metric2(1:end-1).*cos(linspace(-180,179,length(rm2(1,:,1))-1)); % average
        ys = metric2(1:end-1).*sin(linspace(-180,179,length(rm2(1,:,1))-1));
        
        angle2 = atan2(mean(ys),mean(xs)); % mean direction
    end
end
end