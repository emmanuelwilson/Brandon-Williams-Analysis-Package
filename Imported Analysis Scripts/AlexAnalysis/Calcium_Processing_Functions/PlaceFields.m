function [ratemap] = PlaceFields(ms,HD, track, frameMap, thresh)
%%


fps = 30;                                               %Frames per second
spf = 1/fps;                                            %Seconds per frame
ms.timestamp = frameMap.*spf;                           %time stamp in seconds
ms.ind_fire = NaN(ms.numFrames,length(ms.firing(1,:))); %Indices of neuron activity/firing
ms.cell_x = NaN(ms.numFrames,length(ms.firing(1,:)));   %X coordinate of locations where cell fired
ms.cell_y = NaN(ms.numFrames,length(ms.firing(1,:)));   %Y cooridnates of locations where cell fired
ms.cell_time = NaN(ms.numFrames,length(ms.firing(1,:)));%Time at when cell fired
ms.HDfiring = ms.ind_fire;                              %Indices of neuron activity for head direction
mkdir PlaceFieldsRateMap                                    %Create new folder within current directory
direct = dir ('PlaceFieldsRateMap');                         %Access new folder
ratemaps = zeros(37,361,length(ms.firing(1,:)));
probMap = zeros(300,300);
probind = NaN(length(probMap)*length(probMap),1000);

if length(frameMap(:,1)) == length(HD(:,1))
    frameMap2 = zeros(length(frameMap),1);
    for i = 1 : length(frameMap(:,1))
        frameMap2(frameMap(i,1),1) = i;
    end
    frameMap = frameMap2;
end
skipframes = find(~frameMap);

for cellNum = 1 : length(ms.firing(1,:))
    probMap = zeros(300,300);
    probind = NaN(length(probMap)*length(probMap),1000);
    firing = ms.firing(:,cellNum);                          %Extract firing trace
    fire = firing;                                          %Duplicate
    fire(fire < thresh) = 0;                                %apply threshold on firing
    fire(skipframes(:))=0;
    fire(length(track):end) = 0;
    ifire = find(fire);                                     %Find indices for all non-zero values
    
    if(~isempty(ifire))
        for j = 1 : length(ifire)
            ms.ind_fire(j,cellNum) = ifire(j);                      %Add firing index to ms struct
            ms.cell_x(j,cellNum) = track((frameMap(ifire(j))));  %X postion of mouse during firing at sinked time
            ms.cell_y(j,cellNum) = track((frameMap(ifire(j))),2);%Y position of mouse during firing at sinked time
            ms.cell_time(j,cellNum) = ms.timestamp(ifire(j));       %Physiological time of firing
            ms.HDfiring(j,cellNum) = HD(frameMap(ifire(j)));        %Head Direction of mouse at time of neural firing
            if ms.cell_y(j,cellNum) > 0
                n=1;
                while ~isnan(probind(sub2ind(size(probMap),round(ms.cell_y(j,cellNum)-100),round(ms.cell_x(j,cellNum))-100),n))
                    n = n+1;
                end
                probind(sub2ind(size(probMap),round(ms.cell_y(j,cellNum)-100),round(ms.cell_x(j,cellNum))-100),n) = firing(ifire(j));
            end
        end
        probind = mean(probind,2,'omitnan');
        prob = find(probind>0);
        for j = 1 : length(prob)
            probMap(ind2sub(size(probMap),prob(j))) = probind(prob(j));
        end
        
        figure(1);
        n=1;
        c = 1;
        
        %Trajectory map
        subplot(n,2,c);c=c+1;
        hold on       
        plot(track(:,1),-track(:,2),'Color',[.7 .7 .7])
        colormap(gca,hsv)
        xlim([min(track(:,1)) max(track(:,1))]);ylim([-max(track(:,2)) -min(track(:,2))])
        cx=ms.cell_x(:,cellNum);
        cy=ms.cell_y(:,cellNum);
        scatter(cx,-cy,38,ms.HDfiring(:,cellNum),'filled')        
        set(gca,'YDir','Normal')
        title('Traj')
        axis off
        axis square
        hold off
        
        %           Probability heatmap
        subplot(n,2,c);c=c+1;        
        imagesc(probMap)
        colormap(gca,'default')
        colorbar
        
        saveas(gcf,['PlaceFieldsRateMap/',num2str(cellNum),'PlaceField.jpg']); %saving figure as a picture file (.jpg) in the new folder "EBCresults"
        ms.ind_fire = NaN(ms.numFrames,length(ms.firing(1,:))); %Indices of neuron activity/firing
        ms.cell_x = NaN(ms.numFrames,length(ms.firing(1,:)));   %X coordinate of locations where cell fired
        ms.cell_y = NaN(ms.numFrames,length(ms.firing(1,:)));   %Y cooridnates of locations where cell fired
        ms.cell_time = NaN(ms.numFrames,length(ms.firing(1,:)));%Time at when cell fired
        ms.HDfiring = ms.ind_fire;                              %Indices of neuron activity for head direction
        clf
    end
end