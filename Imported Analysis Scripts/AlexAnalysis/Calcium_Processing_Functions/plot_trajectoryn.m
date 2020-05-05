%% Execute this after sync_ms_behav
    a = load('HeadTrackingData.mat');
    b = load('frameMap.mat');
    c = load('trace.mat');
    %d = load('first_peak_time_13_06')
    position_track = a.SINKdata;
    syncedFrames = b.frameMap;
    ms = c.ms;
    %first_peak_time = d.first_peak_time;
for cell_number=1:(length(ms.trace(1,:)))
    
    shifted_trace = [];
    shifted_position_track = [];
    rise_time = 10; % Number of frames on average (try between 20 and 30)
    %cell_number = 1; % number 7, 10 (case study to reduce effect of values below mean fluorescence)
    
    trace = ms.trace(:,cell_number);
    if syncedFrames(length(syncedFrames)) < length(trace)
        trace = trace(syncedFrames);
    end
    firing = ms.firing(:,cell_number);
    
    shifted_trace = trace((rise_time : length(trace)));
    shifted_position_track = (position_track(1 : (length(position_track) - rise_time + 1)))';
    shifted_trace(shifted_trace < 0.1) = 0;
    ind_firing_on_trajectory = find(shifted_trace);
    
    % Uncomment to see firing on trajectory plots
    figure;
    subplot(1, 2, 1)
    plot(position_track(:,1), position_track(:,2))
    hold on
    scatter(position_track(ind_firing_on_trajectory,1), position_track(ind_firing_on_trajectory,2),'r','filled', 'LineWidth', 0.1)
    %polar(sorted_HD, smooth(smooth(sorted_shifted_trace)));
    %polar(sorted_HD_no_shift, (ms.firing(:, cell_number))');
    subplot(1, 2, 2)
    imshow(uint8(ms.frameMax) + 255*uint8(bwperim(ms.segments(:,:,cell_number)))) % replace ms.meanFrame with ms.frameMax
    %% Saving figures
    if ~exist('HD_cells', 'dir') || ~exist('PC_cells', 'dir')
        mkdir('HD_cells');
        mkdir('PC_cells');
    end
    saveas(gcf, ['PC_cells/', num2str(cell_number), 'trajectory_13_06.jpg']);
    %saveas(gcf, ['HD_cells/', num2str(first_peak_time(cell_number)), 'trajectory_13_06.jpg']);
    
    close all
    
end
    