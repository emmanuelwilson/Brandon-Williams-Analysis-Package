%% Execute this after sync_ms_behav
    a = load('5JH2_20173108_of_2_HD_180_deg_total.mat');
    b = load('5JH2_20173108_of_2_frameMap.mat');
    c = load('5JH2_20173108_of_2.mat');
    %d = load('first_peak_time_07_07');
    HD_180_deg = a.HD_180_deg_total';
    syncedFrames = b.frameMap;
    ms = c.ms;
    %first_peak_time = d.first_peak_time;
for cell_number=1:(length(ms.trace(1,:)))
    
    shifted_trace = [];
    shifted_HD_180_deg = [];
    sorted_shifted_trace = [];
    sorted_HD = [];
    ind_sorted_HD = [];
    rise_time = 10; % Number of frames on average (try between 20 and 30)
    %cell_number = 1; % number 7, 10 (case study to reduce effect of values below mean fluorescence)
    
    
%     a = load('HD_180_deg_17_03.mat');
%     b = load('frameMap_17_03.mat');
%     c = load('ms_17_03.mat');
%     HD_180_deg = a.HD_180_deg';
%     syncedFrames = b.frameMap;
%     ms = c.ms;
    
    
    trace = ms.trace(:,cell_number);
    trace = trace(syncedFrames);
    firing = ms.firing(:,cell_number);
    firing = firing(syncedFrames);
    %% Uncomment to see figure
%     figure;plot(trace)
%      hold on
%      plot(HD_180_deg/180)
%      hold on
%      plot(firing)
     
    
    shifted_trace = trace((rise_time : length(trace)));
    shifted_HD_180_deg = (HD_180_deg(1 : (length(HD_180_deg) - rise_time + 1)))';
    
    [sorted_HD, ind_sorted_HD] = sort(deg2rad(shifted_HD_180_deg));
    [sorted_HD_no_shift, ind_sorted_HD_no_shift] = sort(deg2rad(HD_180_deg));
    sorted_shifted_trace = shifted_trace(ind_sorted_HD);
    sorted_shifted_trace(sorted_shifted_trace < 0) = 0;
    firing = firing';      %%%%%%%%%%%%%%%%%%
    firing = firing(ind_sorted_HD_no_shift);    %%%%%%%%%%%%%%%%%%% Synchronisation ??
    sorted_shifted_trace(sorted_shifted_trace < 0) = 0;
    %% Uncomment to see firing polar plots (without binning)
%     figure;
%     subplot(1, 2, 1)
%     polar(sorted_HD_no_shift, firing');
%     subplot(1, 2, 2)
%     imshow(uint8(ms.frameMax) + 255*uint8(bwperim(ms.segments(:,:,cell_number)))) % replace ms.meanFrame with ms.frameMax
%     
    % Uncomment to see calcium polar plots
    figure;
    subplot(1, 2, 1)
    polar(sorted_HD', (smooth(sorted_shifted_trace)));
    %polar(sorted_HD, smooth(smooth(sorted_shifted_trace)));
    %polar(sorted_HD_no_shift, (ms.firing(:, cell_number))');
    subplot(1, 2, 2)
    imshow(uint8(ms.frameMax) + 255*uint8(bwperim(ms.segments(:,:,cell_number)))) % replace ms.meanFrame with ms.frameMax
    
    %% Binning and MRL for sorted shifted traces
    binned_sorted_HD = [];
    binned_sorted_shifted_trace = [];
    sorted_HD = floor(rad2deg(sorted_HD));
    ang_bin_size = 1;
    i = 1;
    temp_sorted_HD = [];
    temp_sorted_shifted_trace = [];
    for ind = 1 : 360/ang_bin_size
    k = 1;
    while (i < (length(sorted_HD) + 1)) && (sorted_HD(i) < (ang_bin_size*ind + sorted_HD(1)))
        temp_sorted_HD(k) = sorted_HD(i);
        temp_sorted_shifted_trace(k) = sorted_shifted_trace(i);
        i = i + 1;
        k = k + 1;
    end
    binned_sorted_HD(ind) = mean(temp_sorted_HD);
    binned_sorted_shifted_trace(ind) = mean(temp_sorted_shifted_trace);
    temp_sorted_HD = [];
    temp_sorted_shifted_trace = [];
    end
    
    binned_sorted_HD = deg2rad(binned_sorted_HD);
    binned_sorted_HD(isnan(binned_sorted_HD)) = [];
    binned_sorted_shifted_trace(isnan(binned_sorted_shifted_trace)) = [];
    figure; 
    subplot(1, 2, 1)
    polar(binned_sorted_HD, ((binned_sorted_shifted_trace)));
    hold on
    
    %% Plotting mean resultant length (need to normalize data over all directions first)
    
    %sorted_HD = sorted_HD * unitsratio('rad','deg');
    
    binned_sorted_HD = cat(1, binned_sorted_HD(:), binned_sorted_HD(1));
    
    binned_sorted_shifted_trace = cat(1, (binned_sorted_shifted_trace(:)), (binned_sorted_shifted_trace(1)));
    %
    % h1=polar(shifted_HD_180_deg(:), shifted_trace(:), 'k'); hold on;
    
    %set(h1,'linewidth',1.1)
    
    xs = binned_sorted_shifted_trace(1:end-1).*cos(binned_sorted_HD(1:end-1)); % average
    ys = binned_sorted_shifted_trace(1:end-1).*sin(binned_sorted_HD(1:end-1));
    
    coordlims=axis;
    
    ang_hd(cell_number) = atan2(mean(ys),mean(xs)); % mean direction
    
    mr = (cos(ang_hd(cell_number))*sum(xs) + sin(ang_hd(cell_number))*sum(ys)) / sum(shifted_trace(1:end-1)); % mean resultant length
    
    mag_hd = sqrt(sum(ys)^2+sum(xs)^2)/sqrt(sum(abs(ys))^2+sum(abs(xs))^2)*coordlims(2); % for visualizations sake
    
    % figure; polar(shifted_HD_180_deg(:), shifted_trace(:), 'k');
    polar([ang_hd(cell_number) ang_hd(cell_number)], [0 mag_hd], 'r')
    subplot(1, 2, 2)
    imshow(uint8(ms.frameMax) + 255*uint8(bwperim(ms.segments(:,:,cell_number)))) % replace ms.meanFrame with ms.frameMax

    
%     %% Binning and MRL for firing
%     binned_sorted_HD_no_shift = [];
%     binned_firing = [];
%     sorted_HD_no_shift = floor(rad2deg(sorted_HD_no_shift));
%     ang_bin_size = 5;
%     i = 1;
%     temp_sorted_HD_no_shift = [];
%     temp_firing = [];
%     for ind = 1 : 360/ang_bin_size
%         k = 1;
%         while (i < (length(sorted_HD_no_shift) + 1)) && (sorted_HD_no_shift(i) < (ang_bin_size*ind + sorted_HD_no_shift(1)))
%             temp_sorted_HD_no_shift(k) = sorted_HD_no_shift(i);
%             temp_firing(k) = firing(i);
%             i = i + 1;
%             k = k + 1;
%         end
%         binned_sorted_HD_no_shift(ind) = mean(temp_sorted_HD_no_shift);
%         binned_firing(ind) = mean(temp_firing);
%         temp_sorted_HD_no_shift = [];
%         temp_firing = [];
%     end
%     
%     binned_sorted_HD_no_shift = deg2rad(binned_sorted_HD_no_shift);
%     binned_sorted_HD_no_shift(isnan(binned_sorted_HD_no_shift)) = [];
%     binned_firing(isnan(binned_firing)) = [];
%     figure;
%     subplot(1, 2, 1)
%     polar(binned_sorted_HD_no_shift, ((binned_firing)));
%     
%     
%     %% Plotting mean resultant length (need to normalize data over all directions first)
%     
%     hold on
%     
%     %sorted_HD = sorted_HD * unitsratio('rad','deg');
%     
%     binned_sorted_HD_no_shift = cat(1, binned_sorted_HD_no_shift(:), binned_sorted_HD_no_shift(1));
%     
%     binned_firing = cat(1, (binned_firing(:)), (binned_firing(1)));
%     %
%     % h1=polar(shifted_HD_180_deg(:), shifted_trace(:), 'k'); hold on;
%     
%     %set(h1,'linewidth',1.1)
%     
%     xs = binned_firing(1:end-1).*cos(binned_sorted_HD_no_shift(1:end-1)); % average
%     ys = binned_firing(1:end-1).*sin(binned_sorted_HD_no_shift(1:end-1));
%     
%     coordlims=axis;
%     
%     ang_hd(cell_number) = atan2(mean(ys),mean(xs)); % mean direction
%     
%     mr = (cos(ang_hd(cell_number))*sum(xs) + sin(ang_hd(cell_number))*sum(ys)) / sum(binned_firing(1:end-1)); % mean resultant length
%     
%     mag_hd = sqrt(sum(ys)^2+sum(xs)^2)/sqrt(sum(abs(ys))^2+sum(abs(xs))^2)*coordlims(2); % for visualizations sake
%     
%     % figure; polar(shifted_HD_180_deg(:), shifted_trace(:), 'k');
%     hold on;
%     polar([ang_hd(cell_number) ang_hd(cell_number)], [0 mag_hd], 'r')
%     subplot(1, 2, 2)
%     imshow(uint8(ms.frameMax) + 255*uint8(bwperim(ms.segments(:,:,cell_number)))) % replace ms.meanFrame with ms.frameMax
    
    %% Saving figures
    if ~exist('HD_cells', 'dir')
        mkdir('HD_cells');
    end
    saveas(gcf, ['HD_cells/', num2str(cell_number), 'Calcium_polar_20_06.jpg']);
    %saveas(gcf, ['HD_cells/', num2str(first_peak_time(cell_number)), 'Calcium_trace_20_06.jpg']);
    
    close all
    
end

% %% Manually choosing most probable HD cells
% most_probable_HD_cells = zeros((length(ms.trace(1,:))), 1);
% most_probable_HD_cells([4 5 10 12 13 15 16 17 18 19 20 21 22 23 24 25 30 32 33 37 38 39 40 42 44 45 46 48 52 53 54 57 58 59 60 64 65 66 67 68 69 70 72 74 76 80 82 84 87 93 95 97 101 107 108 109 110 114 115 116])=1;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Sorting following first peak occurence
% 
% unsorted_first_peaks = [first_peak_time(find(first_peak_time)')' find(first_peak_time)'];
% sorted_first_peaks = sortrows(unsorted_first_peaks);
% 
% for ind1 = 1 : length(sorted_first_peaks)
%     figure;
%     subplot(length(sorted_first_peaks), 1, ind1);
%     plot(ms.trace(:, sorted_first_peaks(ind1,2)));
% end
% 
% for ind2 = 1 : length(sorted_first_peaks)
%     [max_peak(ind2) ind_max_peak(ind2)] = max(ms.trace(:, sorted_first_peaks(ind2,2)));
% end
% %% Heatmap matrix 
% heatmap_matrix = [];
% for ind4 = 1 : length(sorted_first_peaks)
%     heatmap_matrix = [heatmap_matrix ; ms.trace(:, sorted_first_peaks(ind4,2))'];
% end
% figure;
% colormap('hot');
% imagesc(heatmap_matrix);
% colorbar;
% 
% %% Sorting following highest peak occurence
% 
% max_peak = [max_peak' ind_max_peak' sorted_first_peaks(:,2)];
% sorted_max_peak = sortrows(max_peak,2);
% 
% figure;
% for ind3 = 1 : length(sorted_max_peak)
%     subplot(length(sorted_max_peak), 1, ind3);
%     plot(ms.trace(:, sorted_max_peak(ind3,3)));
% end
% %% Heatmap matrix 
% heatmap_matrix = [];
% for ind4 = 1 : length(sorted_max_peak)
%     heatmap_matrix = [heatmap_matrix ; ms.trace(:, sorted_max_peak(ind4,3))'];
% end
% figure;
% colormap('hot');
% imagesc(heatmap_matrix);
% colorbar;


%% Sequence mapping
% fig = uint8(ms.frameMax);
% for cellnum = 1:length(sorted_max_peak(:,3))
%         fig = fig + 255*uint8((ms.segments(:,:,sorted_max_peak(cellnum,3)))); % replace ms.meanFrame with ms.frameMax
% end
% %figure; imshow(fig)
% 
% % Color mapping of HD cells in ROI
% sizeMask = size(ms.mask(:,:));
% 
% %most_probable_HD_cells = ((most_probable_HD_cells.*ang_hd')+pi)/(2*pi);
% %most_probable_HD_cells(most_probable_HD_cells==0.5)=0;
% normalized_time = sorted_max_peak(:,2)/length(ms.trace(:,1));
% map_fig = zeros(sizeMask(1), sizeMask(2),3);
% for cellnum = 1: length(sorted_max_peak(:,3))
%     fig1 = zeros(sizeMask(1), sizeMask(2));
%     hue = zeros(sizeMask(1), sizeMask(2));
%     saturation = zeros(sizeMask(1), sizeMask(2));
%     brightness = zeros(sizeMask(1), sizeMask(2));
%     rgb = zeros(sizeMask(1), sizeMask(2), 3);
%     figfig = zeros(sizeMask(1), sizeMask(2), 3);
%     %if most_probable_HD_cells(cellnum) ~= 0
%         fig1= 255*uint8((ms.segments(:,:,sorted_max_peak(cellnum,3))));
%         fig1 = im2bw(fig1)*normalized_time(cellnum);
%         hue = fig1;
%         saturation = ones(size(hue));
%         brightness = ones(size(hue));
%         rgb = hsv2rgb(cat(3, hue, saturation, brightness));
%         figfig = cat(3, rgb(:,:,1).*im2bw(255*uint8((ms.segments(:,:,sorted_max_peak(cellnum,3))))), rgb(:,:,2).*im2bw(255*uint8((ms.segments(:,:,sorted_max_peak(cellnum,3))))), rgb(:,:,3).*im2bw(255*uint8((ms.segments(:,:,sorted_max_peak(cellnum,3))))));
%     %end
%     map_fig = map_fig + figfig;
% end
% 
% 
% 
% % Plotting Colorbar
% 
% 
% 
% figure
% % positionVector1 = [0.07, 0.15, 0.7, 0.7];
% % subplot('Position',positionVector1)
% imshow(map_fig)




