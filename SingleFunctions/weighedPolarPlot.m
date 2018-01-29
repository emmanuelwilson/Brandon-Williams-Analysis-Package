%% Execute this after sync_ms_behav
a = load('HeadTrackingData.mat');
b = load('frameMap.mat');
c = load('5JH2_20171110_of_1.mat');
syncedFrames = b.frameMap;
HD_180_deg = a.HDdeg(syncedFrames)';
ms = c.ms;
for cell_number=1:(length(ms.trace(1,:)))
    
    shifted_trace = [];
    shifted_HD_180_deg = [];
    sorted_shifted_trace = [];
    sorted_HD = [];
    ind_sorted_HD = [];
    rise_time = 5; % Number of frames on average (try between 20 and 30)
    
    trace = ms.trace(:,cell_number);
    firing = ms.firing(:,cell_number);
    if syncedFrames(length(syncedFrames)) < length(trace)
        trace = trace(syncedFrames);
        firing = firing(syncedFrames);
    end
    
    shifted_trace = trace((rise_time : length(trace)));
    shifted_HD_180_deg = (HD_180_deg(1 : (length(HD_180_deg) - rise_time + 1)))';
    
    [sorted_HD, ind_sorted_HD] = sort(deg2rad(shifted_HD_180_deg));
    [sorted_HD_no_shift, ind_sorted_HD_no_shift] = sort(deg2rad(HD_180_deg));
    sorted_shifted_trace = shifted_trace(ind_sorted_HD);
    sorted_shifted_trace(sorted_shifted_trace < 0) = 0;
    firing = firing';      %%%%%%%%%%%%%%%%%%
    firing = firing(ind_sorted_HD_no_shift);    %%%%%%%%%%%%%%%%%%% Synchronisation ??
    sorted_shifted_trace(sorted_shifted_trace < 0) = 0;
    
    ind1 = 1;
    ind3 = 1;
    ang_bin = 2;
    time_sorted_binned_HD = [];
    aa = [];
    aa = [rad2deg(sorted_HD); ind_sorted_HD]';
    bb = [];
    for ind2 = 1 : 360/ang_bin
        ind3 = 1;
        temp_aa = [];
        sorted_temp_aa = [];
        ind_sort_temp_aa = [];
        while aa(ind1) <= (180 - (360 - ang_bin*ind2)) % Correct this to include cases where animal doesn't span all directions
            temp_aa(ind3, :) = [((ang_bin*ind2-ang_bin)-180) aa(ind1, 2)];
            ind1 = ind1 + 1;
            ind3 = ind3 + 1;
        end
        [sorted_temp_aa(:, 2), ind_sort_temp_aa] = sort(temp_aa(:, 2));
        sorted_temp_aa(:, 1) = temp_aa(ind_sort_temp_aa);
        time_sorted_binned_HD = [time_sorted_binned_HD; sorted_temp_aa];
    end
    
    cc = [];
    time_sorted_binned_shifted_trace = shifted_trace(time_sorted_binned_HD(:, 2));
    bb = [time_sorted_binned_HD time_sorted_binned_shifted_trace];
    cc(1,1) = 1;
    cc(1,2) = cc(1,1) * bb(1,3);
    tolerated_frame_skips = 10;
    for ind4 = 2: length(bb(:,1))
        if (bb((ind4), 2)-bb(ind4-1, 2) < tolerated_frame_skips)
            cc(ind4,1) = cc((ind4-1),1) + 1;
            cc(ind4,2) = cc(ind4,1) * bb(ind4,3);
        else
            cc(ind4,1) = 1;
            cc(ind4,2) = cc(ind4,1) * bb(ind4,3);
        end
    end
    
    dd = [];
    dd = [bb cc];
    
    weighed_polar = [];
    ind5 = 1;
    ind6 = 1;
    ind7 = 1;
    while (ind5 <= length(dd(:,1)))
        ind6 = 1;
        while ((ind5 + ind6)<=length(dd(:,1))) && (dd((ind5),1) == dd((ind5 + ind6) ,1))
            ind6 = ind6 + 1;
        end
        weighed_polar(ind7) = sum(dd((ind5:(ind5+ind6-1)),5))/sum(dd((ind5:(ind5+ind6-1)),4));
        ind5 = ind5 + ind6;
        ind7 = ind7 + 1;
    end
    weighed_polar(weighed_polar < 0) = 0;
    binned_HD = deg2rad([0:ang_bin:359] - 180);
    
    figure;
    subplot(1, 2, 1)
    polar(binned_HD, weighed_polar);
    
    
    %% Plotting mean resultant length (need to normalize data over all directions first)
    
    hold on
    
    %sorted_HD = sorted_HD * unitsratio('rad','deg');
    
    binned_HD = cat(1, binned_HD(:), binned_HD(1));
    
    weighed_polar(weighed_polar < 0.02) = 0; % change conditions accordingly (ex: < 0.05)
    weighed_polar = cat(1, (weighed_polar(:)), (weighed_polar(1)));
    %
    % h1=polar(shifted_HD_180_deg(:), shifted_trace(:), 'k'); hold on;
    
    %set(h1,'linewidth',1.1)
    
    xs = weighed_polar(1:end-1).*cos(binned_HD(1:end-1)); % average
    ys = weighed_polar(1:end-1).*sin(binned_HD(1:end-1));
    
    coordlims=axis;
    
    ang_hd(cell_number) = atan2(mean(ys),mean(xs)); % mean direction
    
    mr = (cos(ang_hd(cell_number))*sum(xs) + sin(ang_hd(cell_number))*sum(ys)) / sum(weighed_polar(1:end-1)); % mean resultant length
    
    mag_hd = sqrt(sum(ys)^2+sum(xs)^2)/sqrt(sum(abs(ys))^2+sum(abs(xs))^2)*coordlims(2); % for visualizations sake
    
    % figure; polar(shifted_HD_180_deg(:), shifted_trace(:), 'k');
    polar([ang_hd(cell_number) ang_hd(cell_number)], [0 mag_hd], 'r')
    subplot(1, 2, 2)
    imshow(uint8(ms.frameMax) + 255*uint8(bwperim(ms.segments(:,:,cell_number)))) % replace ms.meanFrame with ms.frameMax     
    
    
    %% Saving figures
    if ~exist('HD_cells', 'dir')
        mkdir('HD_cells');
    end
    saveas(gcf, ['HD_cells/', num2str(cell_number), '_weighed_polar_ab3_rt5_12705_11_07.jpg']); % ab: angle bin / rt: rise time
    
    close all
    
end


% %% Manually choosing most probable HD cells
% most_probable_HD_cells = zeros((length(ms.trace(1,:))), 1);
% most_probable_HD_cells([1 2 3 6 16 17 22 25 26 28 39 42 43 47 52 59 61 67 69 74 76 77 79 80 86 88 89 92 93 94 103 105 110 115])=1;
% fig = uint8(ms.frameMax);
% for cellnum = 1:length(most_probable_HD_cells)
%     if most_probable_HD_cells(cellnum) == 1
%         fig = fig + 255*uint8((ms.segments(:,:,cellnum))); % replace ms.meanFrame with ms.frameMax
%     end
% end
% %figure; imshow(fig)
% 
% % Color mapping of HD cells in ROI
% sizeMask = size(ms.mask(:,:));
% most_probable_HD_cells = ((most_probable_HD_cells.*ang_hd')+pi)/(2*pi);
% most_probable_HD_cells(most_probable_HD_cells==0.5)=0;
% map_fig = zeros(sizeMask(1), sizeMask(2),3);
% for cellnum = 1: length(most_probable_HD_cells)
%     fig1 = zeros(sizeMask(1), sizeMask(2));
%     hue = zeros(sizeMask(1), sizeMask(2));
%     saturation = zeros(sizeMask(1), sizeMask(2));
%     rgb = zeros(sizeMask(1), sizeMask(2), 3);
%     figfig = zeros(sizeMask(1), sizeMask(2), 3);
%     if most_probable_HD_cells(cellnum) ~= 0
%         fig1= 255*uint8((ms.segments(:,:,cellnum)));
%         fig1 = im2bw(fig1)*most_probable_HD_cells(cellnum);
%         hue = fig1;
%         saturation = ones(size(hue));
%         brightness = ones(size(hue));
%         rgb = hsv2rgb(cat(3, hue, saturation, brightness));
%         figfig = cat(3, rgb(:,:,1).*im2bw(255*uint8((ms.segments(:,:,cellnum)))), rgb(:,:,2).*im2bw(255*uint8((ms.segments(:,:,cellnum)))), rgb(:,:,3).*im2bw(255*uint8((ms.segments(:,:,cellnum)))));
%     end
%     map_fig = map_fig + figfig;
% end
% 
% 
% 
% % Plotting Color key
% % Set parameters (these could be arguments to a function)
% rInner = 80;     % inner radius of the colour ring
% rOuter = 200;    % outer radius of the colour ring
% % Get polar coordinates of each point in the domain
% [x, y] = meshgrid(-rOuter:rOuter);
% [theta, rho] = cart2pol(x, -y);
% % Set up colour wheel in hsv space
% hue = (theta + pi) / (2 * pi);     % hue into range (0, 1]
% saturation = ones(size(hue));      % full saturation
% brightness = double(rho >= rInner & rho <= rOuter);  % black outside ring
% % Convert to rgb space for display
% rgb = hsv2rgb(cat(3, hue, saturation, brightness));
% 
% 
% 
% figure
% positionVector1 = [0.07, 0.15, 0.7, 0.7];
% subplot('Position',positionVector1)
% imshow(map_fig)
% 
% positionVector2 = [0.77, 0.4, 0.2, 0.2];
% subplot('Position',positionVector2)
% imshow(rgb)



%% Show frameMax
% figure;imshow(imadjust(ms.frameMax, [0 0.2], [0 1]))