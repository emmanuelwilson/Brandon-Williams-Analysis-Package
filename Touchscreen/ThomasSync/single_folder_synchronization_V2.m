function single_folder_synchronization_V2(data_dir)

% Function that reads the timestamp.csv, schedules.csv files and videos
% labeled msCam[].avi and behaviorCam[].avi, and builds synchronization
% matrices that are saved in msTouchSynch.mat.

% timestamp.csv: file produced by the experiment containing times when
% cameras 0 and 1 (measured by the system clock) took frames, which were
% then saved in video files.

% schedules.csv: file that contains the output from the cage during the
% recording session, with time measured relative to the cage clock.

% msTouchSynch.mat: contains the structure synchronization with the
% following substructures:
% synchronization.miniscopeMaster = synchronization matrix with miniscope
% time as reference time. Frames and times for miniscope and behavior are
% held in synchronization.miniscopeMaster.masterTimes/.masterFrames/
% .slaveTimes/.slaveFrames.
% synchronization.behaviorMaster = same but with behavior taken as
% reference.

% 1 = show progress, 0 = do not show
show_progress = 0;
oldcd = pwd;
cd(data_dir)
% timestamp file name
timestamp_file = 'timestamp.dat';
% schedule.csv file
csvFiles = dir('*.csv');
schedule_file = csvFiles(1).name;

nb_lines = 0;

% array to store content of timestamp.dat
camNum_array = [];
frameNum_array = [];
sysClock_array = [];
buffer_array = [];

% number of cam0 and cam1 entries
Ncam0 = 0;
Ncam1 = 0;

% total number of frames within videos
mini_numFrames = 0;
behav_numFrames = 0;
% total number of each video files
mini_numVideos = 0;
behav_numVideos = 0;
% mean signal in each video type
mean_msSignal = [];
mean_bhSignal = [];

% labels of cams 0 and 1 (= miniscope or behavior)
label0 = '';
label1 = '';

% define structure holding content of timestamp.dat
timestamp = [];

% output structure
synchronization = [];

% open figures for later
screensize = get(0,'Screensize');

% window showing progress in reading videos
h_progress = figure('Position',[0.25*screensize(3) screensize(4)/3 screensize(3)/2 0.9*screensize(4)/2],'Visible','Off');

% content of timestamp
h_timestamp = figure('Position',[1 screensize(4)/2 screensize(3) 0.9*screensize(4)/2],'Visible','Off');

% content of videos
h_videos = figure('Position',[1 screensize(4)*1/3 screensize(3) 0.9*screensize(4)*2/3],'Visible','Off');

% comparison of the the two and the synchronization
h_synchr = figure('Position',[1 screensize(4)*1/4 screensize(3) 0.9*screensize(4)*3/4],'Visible','Off');

% flag to skip datasets if there are problems
dataset_ok = 1;

% size of window used to compute the mean frame rate
window_size = 100;

% -------------------------------------------------------------------------



% read the content of timestamp.csv and store it in arrays
[camNum_array,frameNum_array,sysClock_array,buffer_array,Ncam0,Ncam1] = read_timestamp([data_dir '\' timestamp_file]);


% text summary of the raw content of the frames/sysClock values for both
% cameras
if dataset_ok ==1
    disp_timeStpFile(camNum_array,frameNum_array,sysClock_array,Ncam0,Ncam1);
end

% plot system clocks
if dataset_ok ==1
    plot_clocks(camNum_array,sysClock_array,Ncam0,Ncam1);
end

% plot buffers
if dataset_ok ==1
    plot_buffers(camNum_array,buffer_array,Ncam0,Ncam1);
end

% buffer content inspection - check for large number of frames kept in
% buffer that might indicate discarded frames
if dataset_ok ==1
    disp_buffer(camNum_array,buffer_array)
end

% SysClock might be equal to large negative value(s) for the initial
% frame(s). One is ok and might signify a clock value almost equal to 0.
% More should be investigated.
if dataset_ok ==1
    [sysClock_array] = remove_faulty_systemclock(camNum_array,sysClock_array);
end

% plot also the frame rate for both cameras
if dataset_ok ==1
    disp_frame_rate(camNum_array,frameNum_array,sysClock_array,buffer_array);
end

if dataset_ok ==1
    % count the number of frames in the miniscope and behavioral videos
    [mini_numVideos,mini_numFrames,behav_numVideos,behav_numFrames,mean_msSignal,mean_bhSignal] = read_videos(data_dir,Ncam0,Ncam1);
    
    % display the number of frames in each video and the mean intensity
    disp_video_stats(mini_numVideos,mini_numFrames,behav_numVideos,behav_numFrames,mean_msSignal,mean_bhSignal);
end

% display number of frames in cams 0, 1, miniscope and behavior and how
% they are matched
if dataset_ok ==1
    display_nbFrames(Ncam0,Ncam1,mini_numFrames,behav_numFrames,mini_numVideos,behav_numVideos);
end

% fill timestamp structure containing all the relevant information so far
if dataset_ok ==1
    [timestamp] = fill_timestamp(camNum_array,frameNum_array,sysClock_array,buffer_array,timestamp);
end

% plot of the frames from both cameras as as function of time and a linear
% fit to see which is reliable implementation of time in the experiment
if dataset_ok ==1
    disp_time_vs_frames(timestamp,mini_numFrames,behav_numFrames)
end

% synchronize the three streams of information
% there are two possibilities:
% 1) synchronizing the behavior to the miniscope (miniscope = master, behavior = slave)
synchronization.miniscopeMaster.masterFrames = [];
synchronization.miniscopeMaster.masterTimes = [];
synchronization.miniscopeMaster.slaveFrames = [];
synchronization.miniscopeMaster.slaveTimes = [];
synchronization.miniscopeMaster.N = [];

% 2) synchronizing the miniscope to the behavior (behavior = master, miniscope = slave).
synchronization.behaviorMaster.masterFrames = [];
synchronization.behaviorMaster.masterTimes = [];
synchronization.behaviorMaster.slaveFrames = [];
synchronization.behaviorMaster.slaveTimes = [];
synchronization.behaviorMaster.N = [];

if dataset_ok ==1
    [synchronization] = synchronize(timestamp);
end

% adding the synchronized cage time
if dataset_ok ==1
    [synchronization] = compute_cage_time(synchronization,[data_dir '\' schedule_file]);
end

% remove redundant entries at the end of the synchronization matrices
if dataset_ok ==1
    [synchronization] = cleanup(synchronization);
end

% display the content of the synchronization matrices
if dataset_ok ==1
    disp_synch_matrices(synchronization,timestamp);
end

% save figures as diagnostics
save_figures;

% Add a substructure event that contains the events object that Emmanuel's
% code needs to sort out delay cell activity. This makes the assumption
% that the miniscope is the master time series.

% The 6 columns of events are as follows:
% events(:,1) = frameMap (behavior frames)
% events(:,2) = timestamp of event (behavior times)
% events(:,3) = event name
% events(:,4) = nose-poke position
% events(:,5) = Trial Start/stop, indicated by a 1
% events(:,6) = Delay period start/stop (1 = start, 2 = stop)
if dataset_ok ==1
    [synchronization] = add_events(synchronization,[data_dir '\' schedule_file]);
end

% finally, store the results
if dataset_ok ==1
    disp(['Saving the synchronization matrix in file ' data_dir '\msTouchSync_new.mat']);
    save([data_dir '\msTouchSync_new.mat'],'synchronization');
    disp('Done');
end

cd(oldcd)
% ================================= functions =============================

% read the content of timestamp.csv and store it in arrays
    function [cam,frame,cl,buff,N0,N1] = read_timestamp(fname)
        
        cam = [];
        frame = [];
        cl = [];
        buff = [];
        
        % read file timestamp.csv -------------------------------------
        fid = fopen(fname,'r');         %access dat files
        dataArray = textscan(fid, '%f%f%f%f%[^\n\r]', 'Delimiter', '\t', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);    %read file and make sure it is not empty
        cam = dataArray{:, 1};       %camera number
        frame = dataArray{:, 2};     %frame number
        cl = dataArray{:, 3};     %system clock
        buff = dataArray{:, 4};      %buffer
        clearvars dataArray;            %clear variables from dataArray
        fclose(fid);
        
        % check its content ---------------------------------------------
        % List the number of cameras. There should be only two, one labeled 0 and
        % the other 1.
        if ~(min(unique(cam)) == 0 && max(unique(cam)) == 1)
            disp('File timestamp.dat does not feature the required cameras 0 and 1. Stopping here...');
            dataset_ok = 0;
            
        end
        
        % store number of frames for each camera
        N0 = length(find(cam==0));
        N1 = length(find(cam==1));
    end

% text summary of the raw content of the frames/sysClock values for both
% cameras
    function disp_timeStpFile(cam,frame,cl,N0,N1)
        
        f = figure(h_timestamp);
        f.Visible = 'On';
        subplot(1,6,1)
        temp0 = find(cam==0);
        tframes0 = frame(temp0);
        ttimes0 = cl(temp0);
        
        temp1 = find(cam==1);
        tframes1 = frame(temp1);
        ttimes1 = cl(temp1);
        
        % display stuff using text()
        fsize = 8;
        dx = 2.5*2;
        dy = 1*2;
        
        nb = 20;
        
        axis([-dx/2 4*dx -(nb+2)*dy dy])
        
        % column labels
        y = 0;
        text(0,y,'frames 0','FontSize',fsize)
        text(1*dx,y,'times 0','FontSize',fsize)
        text(2*dx,y,'frames 1','FontSize',fsize)
        text(3*dx,y,'times 1','FontSize',fsize)
        
        % matrix content
        for i1=1:nb
            y = y - dy;
            text(0,y,num2str(tframes0(i1)),'FontSize',fsize)
            text(1*dx,y,num2str(ttimes0(i1)),'FontSize',fsize)
            text(2*dx,y,num2str(tframes1(i1)),'FontSize',fsize)
            text(3*dx,y,num2str(ttimes1(i1)),'FontSize',fsize)
        end
        
        % ...
        y = y - dy;
        text(0,y,'...','FontSize',fsize)
        text(1*dx,y,'...','FontSize',fsize)
        text(2*dx,y,'...','FontSize',fsize)
        text(3*dx,y,'...','FontSize',fsize)
        
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        
        xlabel(['Camera 0: ' num2str(N0) ' frames, Camera 1: ' num2str(N1) ' frames.']);
        
        drawnow
        
    end


% displays sysClock for both cameras
    function plot_clocks(cam,cl,N0,N1)
        
        f = figure(h_timestamp);
        subplot(1,6,2)
        
        % single out the outliers
        temp0 = cl(cam==0);
        temp1 = cl(cam==1);
        maxend = max([temp0(end) temp1(end)]);
        
        hold on
        plot(1:length(temp0),cl(cam==0),'m-')
        plot(1:length(temp1),cl(cam==1),'c-')
        axis([-0.1*max([N0 N1]) max([N0 N1]) -1.1*maxend 1.1*maxend])
        grid on
        ylabel(['system clock, max(0) = ' num2str(temp0(end)/1000) 's, max(1) = ' num2str(temp1(end)/1000) 's'])
        xlabel('index, cam0 (mauve)/cam1 (cyan)')
    end


% display buffer for both cameras
    function plot_buffers(cam,buff,N0,N1)
        
        f = figure(h_timestamp);
        subplot(1,6,3)
        title(['Data Processed is in folder ' data_dir],'Interpreter','None');
        
        % single out the outliers
        temp0 = buff(cam==0);
        temp1 = buff(cam==1);
        
        hold on
        plot(temp0,'m-')
        plot(temp1,'c-')
        mmax = max([length(temp0) length(temp1)]);
        axis([0 mmax 0 10])
        grid on
        ylabel('buffer')
        xlabel('index, cam0 (mauve)/cam1 (cyan)')
    end


% buffer content inspection - check for large number of frames kept in
% buffer that might indicate discarded frames
    function disp_buffer(cam,buff)
        
        % We plot the buffer values separately for cam0 and cam1
        % cam0
        temp = find(cam==0);
        all_values = unique(buff(temp));
        
        f = figure(h_timestamp);
        subplot(1,6,4)
        hold on
        labels = [];
        for i1=1:length(all_values)
            my_val = all_values(i1);
            if my_val>1
                bar(i1,length(find(buff(temp)==my_val)),'m')
                labels{end+1} = num2str(my_val);
            end
        end
        set(gca,'XTick',2:length(all_values),'XTickLabel',labels);
        xlabel(['Buffer values (>1). Buffer - Frames = ' num2str(sum(buff) - size(cam,2))]);
        ylabel('Frequencies for cam 0')
        
        % cam1
        temp = find(cam==1);
        all_values = unique(buff(temp));
        
        subplot(1,6,5)
        hold on
        labels = [];
        for i1=1:length(all_values)
            my_val = all_values(i1);
            if my_val>1
                bar(i1,length(find(buff(temp)==my_val)),'c')
                labels{end+1} = num2str(my_val);
            end
        end
        set(gca,'XTick',2:length(all_values),'XTickLabel',labels);
        xlabel(['Buffer values (>1). Buffer - Frames = ' num2str(sum(buff) - size(cam,2))]);
        ylabel('Frequencies for cam 1')
        
        drawnow
        
    end


% remove negative system clock entry
    function [newcl] = remove_faulty_systemclock(cam,cl)
        
        % It seems that the first value of the clock for one camera or the
        % other is large (though sometimes not large enough to be the
        % maximum of all system clock values). To clean this up, we will
        % look at the first clock value for each camera, and set it to 0 if
        % it is larger than the following value. It can also be negative,
        % so we will use the absolute value.
        cam0 = find(cam==0);
        if abs(cl(cam0(1)))>abs(cl(cam0(2)))
            disp(['First entry of cam 0 is equal to ' num2str(cl(cam0(1))) ' and is large: setting it to 0.']);
            cl(cam0(1)) = 0;
        end
        
        cam1 = find(cam==1);
        if abs(cl(cam1(1)))>abs(cl(cam1(2)))
            disp(['First entry of cam 1 is equal to ' num2str(cl(cam1(1))) ' and is large: setting it to 0.']);
            cl(cam1(1)) = 0;
        end
        
        % output cleaned clock entries
        newcl = cl;
    end


% compute and display the frame rate
    function disp_frame_rate(cam,frame,cl,buff)
        
        % frame rate computed per ... frames
        f = figure(h_timestamp);
        subplot(1,6,6)
        
        % cam 0
        frames0 = frame(find(cam==0));
        N0 = ceil(frames0(end)/window_size);
        times0 = cl(cam==0);
        rate0 = zeros(1,N0);
        
        for i1=1:N0
            if i1<N0
                t_values = times0(window_size*(i1-1)+(1:window_size));
                dt = t_values(end) - t_values(1);
                f_values = frames0(window_size*(i1-1)+(1:window_size));
                df = f_values(end) - f_values(1);
            else
                t_values = times0(window_size*(i1-1):end);
                dt = t_values(end) - t_values(1);
                f_values = frames0(window_size*(i1-1):end);
                df = f_values(end) - f_values(1);
            end
            rate0(i1) = dt/df;
        end
        plot(1:N0,rate0,'-mo')
        % plot buffer on top
        hold on
        b0 = buff(cam==0);
        plot(linspace(1,N0,length(b0)),b0,'Color',[1 0 0])
        
        % cam 1
        frames1 = frame(find(cam==1));
        N1 = ceil(frames1(end)/window_size);
        times1 = cl(cam==1);
        rate1 = zeros(1,N1);
        for i1=1:N1
            if i1<N1
                t_values = times1(window_size*(i1-1)+(1:window_size));
                dt = t_values(end) - t_values(1);
                f_values = frames1(window_size*(i1-1)+(1:window_size));
                df = f_values(end) - f_values(1);
            else
                t_values = times1(window_size*(i1-1):end);
                dt = t_values(end) - t_values(1);
                f_values = frames1(window_size*(i1-1):end);
                df = f_values(end) - f_values(1);
            end
            rate1(i1) = dt/df;
        end
        plot(1:N1,rate1,'-co')
        % plot buffer on top
        hold on
        b1 = buff(cam==1);
        plot(linspace(1,N1,length(b1)),b1,'Color',[0 1 0])
        
        ylabel('1/instantaneous frame rate (s)')
        xlabel('cam 0 (mauve)/cam 1 (cyan)/buffer (red/green)');
        
        temp = union(rate0,rate1);
        ylim([min(temp) max(temp)]);
        
    end


% count the number of frames in the miniscope and behavioral videos
    function [msNvideos,msNframes,bhNvideos,bhNframes,mean_ms,mean_bh] = read_videos(vd,N0,N1)
        
        % find avi and dat files
        aviFiles = dir([vd '\*.avi']);
        
        % extract the numbers of msCam and behavCam files
        msNums = [];
        bhNums = [];
        
        for i1=1:length(aviFiles)
            my_name = aviFiles(i1).name;
            
            % separate miniscope and behavioral video files
            if contains(my_name,'msCam')
                % strip the name to extract the video number
                my_name = my_name(1:end-4);
                my_name = my_name(6:end);
                % and store the number
                msNums = [msNums str2num(my_name)];
                
            elseif contains(my_name,'behavCam')
                % strip the name to extract the video number
                my_name = my_name(1:end-4);
                my_name = my_name(9:end);
                % and store the number
                bhNums = [bhNums str2num(my_name)];
            end
        end
        
        % reorder the numbers
        msNums = unique(msNums);
        if isequal(msNums,1:msNums(end))==0
            disp('Problem with miniscope video numbering. Stopping here.')
            dataset_ok = 0;
            
        end
        
        bhNums = unique(bhNums);
        if isequal(bhNums,1:bhNums(end))==0
            disp('Problem with behavioral video numbering. Stopping here.')
            dataset_ok = 0;
            
        end
        
        % miniscope files
        % ---------------
        msNframes = zeros(1,length(msNums));    % number of frames within miniscope videos
        msNvideos = length(msNums);             % nb of videos
        mean_ms = [];
        
        % display progress on screen
        if show_progress
            f = figure(h_progress);
            f.Visible = 'On';
            subplot(1,2,1)
            hold on
            axis([-1 1 -1 1])
            axis square
            patch(cos(0:0.01:2*pi),sin(0:0.01:2*pi),[1 1 1])
            set(gca,'XTick',[],'YTick',[]);
            xlabel('Loading Miniscope Videos...');
        end
        
        for i1=1:msNvideos
            
            % time how long it takes to load one file
            if i1==1
                tic
            end
            
            avi_name = [vd '\msCam' num2str(msNums(i1)) '.avi'];
            mat_name = [vd '\msCam' num2str(msNums(i1)) '.mat'];
            mean_file = [vd '\ms_mean' num2str(msNums(i1)) '.mat'];
            
            % extract the number of frames in the video
            msvidObj = VideoReader(avi_name);
            msNframes(i1) = msvidObj.NumberOfFrames;
            
            % display the mean intensity in the video. First, look for
            % files ms_mean.mat. If it exists, load it and display
            % its content. If it is not there, load the msCam mat files,
            % which are quicker to load than the avi files. If these are
            % not present either, then load the avi files.
            if isfile(mean_file)
                % mean intensity file
                disp(['Loading file ' mean_file]);
                load(mean_file);
                mean_ms = [mean_ms ms_mean_signal];
                
            elseif isfile(mat_name)
                % mat video file
                disp(['Loading file ' mat_name]);
                load(mat_name);
                % compute the mean intensity
                temp = squeeze(video);
                temp = mean(temp,1);
                temp = squeeze(mean(temp,2))';
                mean_ms = [mean_ms temp];
                % save it
                disp(['Saving file ' mean_file]);
                ms_mean_signal = temp;
                save(mean_file,'ms_mean_signal');
                
            elseif isfile(avi_name)
                % avi video file
                disp(['Loading file ' avi_name]);
                
                try
                    video=msvidObj.read();
                    
                    % compute the number of frames of the video and the mean intensity
                    temp = squeeze(video);
                    temp = mean(temp,1);
                    temp = squeeze(mean(temp,2))';
                    mean_ms = [mean_ms temp];
                    % save it
                    disp(['Saving file ' mean_file]);
                    ms_mean_signal = temp;
                    save(mean_file,'ms_mean_signal');
                    disp('Done');
                catch
                    disp(['Video ' avi_name ' could not be read. Moving on...']);
                end

            end
            
            disp(' ')
            
            % time how long it takes to load one file
            if i1==1
                load_time = toc;
            end
            
            % draw quadrant
            if show_progress
                f = figure(h_progress);
                p = patch([0 cos((2*pi*(i1-1)/msNvideos)+(0:0.01:(2*pi/msNvideos+0.01))) 0],[0 sin((2*pi*(i1-1)/msNvideos)+(0:0.01:(2*pi/msNvideos+0.01))) 0],'r');
                p.EdgeColor = [1 0 0];
                
                tot_time = msNvideos*load_time;
                remain_time = (tot_time-i1*load_time)/60;
                remain_time = round(remain_time*10)/10;
                xlabel(['Loading Miniscope Videos - remaining time is about ' num2str(remain_time) ' minutes']);
                
                drawnow
            end
        end
        
        disp(['Miniscope videos: nb of videos = ' num2str(msNvideos) ' and ' num2str(sum(msNframes)) ' frames in all.'])
        disp(' ')
        
        % behavior files
        % ---------------
        bhNframes = zeros(1,length(bhNums));    % number of frames within miniscope videos
        bhNvideos = length(bhNums);             % nb of videos
        mean_bh = [];
        
        % display progress on screen
        if show_progress
            f = figure(h_progress);
            subplot(1,2,2)
            hold on
            axis([-1 1 -1 1])
            axis square
            patch(cos(0:0.01:2*pi),sin(0:0.01:2*pi),[1 1 1])
            set(gca,'XTick',[],'YTick',[]);
            xlabel('Loading Behavior Videos...');
        end
        
        for i1=1:bhNvideos
            
            % time how long it takes to load one file
            if i1==1
                tic
            end
            
            avi_name = [vd '\behavCam' num2str(bhNums(i1)) '.avi'];
            mat_name = [vd '\behavCam' num2str(bhNums(i1)) '.mat'];
            mean_file = [vd '\bh_mean' num2str(bhNums(i1)) '.mat'];
            
            % extract the number of frames in the video
            bhvidObj = VideoReader(avi_name);
            bhNframes(i1) = bhvidObj.NumberOfFrames;
            
            % display the mean intensity in the video. First, look for
            % files bh_mean.mat. If it exists, load it and display
            % its content. If it is not there, load the bhCam mat files,
            % which are quicker to load than the avi files. If these are
            % not present either, then load the avi files.
            if isfile(mean_file)
                % mean intensity file
                disp(['Loading file ' mean_file]);
                load(mean_file);
                mean_bh = [mean_bh bh_mean_signal];
                
            elseif isfile(mat_name)
                % mat video file
                disp(['Loading file ' mat_name]);
                load(mat_name);
                % compute the mean intensity
                temp = squeeze(video);
                temp = mean(temp,1);
                temp = mean(temp,2);
                temp = squeeze(mean(temp,3))';
                mean_bh = [mean_bh temp];
                % save it
                disp(['Saving file ' mean_file]);
                bh_mean_signal = temp;
                save(mean_file,'bh_mean_signal');
                
            elseif isfile(avi_name)
                % avi video file
                disp(['Loading file ' avi_name]);
                
                try
                    video=bhvidObj.read();
                    
                    % compute the number of frames of the video and the mean intensity
                    temp = squeeze(video);
                    temp = mean(temp,1);
                    temp = mean(temp,2);
                    temp = squeeze(mean(temp,3))';
                    mean_bh = [mean_bh temp];
                    % save it
                    disp(['Saving file ' mean_file]);
                    bh_mean_signal = temp;
                    save(mean_file,'bh_mean_signal');
                    disp('Done');
                catch
                    disp(['Video ' avi_name ' could not be read. Moving on...']);
                end
            end
            
            disp(' ')
            
            % time how long it takes to load one file
            if i1==1
                load_time = toc;
            end
            
            % draw quadrant
            if show_progress
                f = figure(h_progress);
                p = patch([0 cos((2*pi*(i1-1)/bhNvideos)+(0:0.01:(2*pi/bhNvideos+0.01))) 0],[0 sin((2*pi*(i1-1)/bhNvideos)+(0:0.01:(2*pi/bhNvideos+0.01))) 0],'g');
                p.EdgeColor = [0 1 0];
                
                tot_time = bhNvideos*load_time;
                remain_time = (tot_time-i1*load_time)/60;
                remain_time = round(remain_time*10)/10;
                xlabel(['Loading Behavior Videos - remaining time is about ' num2str(remain_time) ' minutes']);
                
                drawnow
            end
        end
        
        disp(['Behavior videos: nb of videos = ' num2str(bhNvideos) ' and ' num2str(sum(bhNframes)) ' frames in all.'])
        disp(' ')
        
        % try to match cams 0 and 1, and that of the miniscope and behavior
        if isequal([N0 N1],[sum(msNframes) sum(bhNframes)])
            label0 = 'miniscope';
            label1 = 'behavior';
            
            disp('Frame number match has been found: cam 0 = miniscope, cam 1 = behavior');
        else
            if isequal([N0 N1],[sum(bhNframes) sum(msNframes)])
                label0 = 'behavior';
                label1 = 'miniscope';
                
                disp('Frame number match has been found: cam 0 = behavior, cam 1 = miniscope');
            else
                disp('The number of frames in cameras 0 and 1 do not match those of cameras miniscope and behavior. Stopping here...')
                
                dataset_ok = 0;
            end
        end
        
        close(h_progress);
        
    end


% display the number of frames in each video, and the mean intensity of
% each
    function disp_video_stats(msNvid,msNfr,bhNvid,bhNfr,mean_ms,mean_bh)
        
        f = figure(h_videos);
        f.Visible = 'On';
        
        subplot(2,5,1)
        
        plot(msNfr,'r-o')
        ylabel('Number of frames');
        xlabel('Miniscope video number');
        axis tight
        
        subplot(2,5,2:5)
        plot(mean_ms,'r')
        ylabel('Mean signal');
        xlabel('Miniscope video frames (by trial) - time derivative (blue)');
        axis tight
        set(gca,'XTick',cumsum(msNfr),'XTickLabel',0:msNvid);
        grid on
        
        title(['Data Processed is in folder ' data_dir],'Interpreter','None');
        
        hold on
        plot(3:length(mean_ms),diff(mean_ms(2:end)))
        
        disp(['Miniscope videos: nb of videos = ' num2str(msNvid) ' and ' num2str(sum(msNfr)) ' frames in all.'])
        
        subplot(2,5,6)
        plot(bhNfr,'g-o')
        ylabel('Number of frames');
        xlabel('Behavior video number');
        axis tight
        
        
        subplot(2,5,7:10)
        plot(mean_bh,'g')
        ylabel('Mean signal');
        xlabel('Behavior video frames (by trial) - time derivative (blue)');
        axis tight
        set(gca,'XTick',cumsum(bhNfr),'XTickLabel',0:bhNvid);
        grid on
        
        hold on
        plot(3:length(mean_bh),diff(mean_bh(2:end)))
        
        disp(['Behavior videos: nb of videos = ' num2str(bhNvid) ' and ' num2str(sum(bhNfr)) ' frames in all.'])
        
        drawnow
        
    end


% display number of frames in cams 0, 1, miniscope and behavior
    function display_nbFrames(N0,N1,msNframes,bhNframes,msNvideos,bhNvideos)
        
        f = figure(h_synchr);
        f.Visible = 'On';
        subplot(2,4,1)
        
        % cams are in the order 0, 1, miniscope and behavior
        if isequal([N0 N1],[sum(msNframes) sum(bhNframes)])
            cam_colors = ['r';'g';'r';'g'];
        elseif isequal([N1 N0],[sum(msNframes) sum(bhNframes)])
            cam_colors = ['g';'r';'r';'g'];
        end
        
        hold on
        
        bar(1,N0,cam_colors(1))
        text(1,1.1*N0,['N = ' num2str(N0)],'HorizontalAlignment','Center');
        
        bar(2,N1,cam_colors(2))
        text(2,1.1*N1,['N = ' num2str(N1)],'HorizontalAlignment','Center');
        
        smsNframes = sum(msNframes);
        sbhNframes = sum(bhNframes);
        
        bar(3,smsNframes,cam_colors(3));
        text(3,1.1*smsNframes,['N = ' num2str(smsNframes)],'HorizontalAlignment','Center');
        h = text(3,0.5*smsNframes,[num2str(msNvideos) ' video files'],'HorizontalAlignment','Center');
        set(h,'Rotation',90);
        
        bar(4,sbhNframes,cam_colors(4));
        text(4,1.1*sbhNframes,['N = ' num2str(sbhNframes)],'HorizontalAlignment','Center');
        h = text(4,0.5*sbhNframes,[num2str(bhNvideos) ' video files'],'HorizontalAlignment','Center');
        set(h,'Rotation',90);
        
        xlabel(['camera 0 = ' label0 ' and camera 1 = ' label1 '.']);
        ylabel('Number of frames')
        xlim([0 5])
        set(gca,'XTick',1:4,'XTickLabel',{'cam0','cam1','miniscope','behavior'});
        
    end


% define structure holding content of timestamp.dat
    function [newts] = fill_timestamp(cam,frame,cl,buff,ts)
        
        % store number of frames for each camera
        times0 = cl(cam==0);
        frames0 = frame(cam==0);
        buffer0 = buff(cam==0);
        times1 = cl(cam==1);
        frames1 = frame(cam==1);
        buffer1 = buff(cam==1);
        
        ts.cam0.N = length(times0);
        ts.cam1.N = length(times1);
        
        % store the frames, buffer values and times for each camera
        ts.cam0.frames = frames0;
        ts.cam0.times = times0;
        ts.cam0.buffer = buffer0;
        ts.cam1.frames = frames1;
        ts.cam1.times = times1;
        ts.cam1.buffer = buffer1;
        
        % camera identities
        ts.cam0.label = label0;
        ts.cam1.label = label1;
        
        % and reverse labels
        if strcmp(ts.cam0.label,'miniscope') && strcmp(ts.cam1.label,'behavior')
            ts.miniscope.label = '0';
            ts.behavior.label = '1';
        else
            if strcmp(ts.cam0.label,'behavior') && strcmp(ts.cam1.label,'miniscope')
                ts.miniscope.label = '1';
                ts.behavior.label = '0';
            else
                disp('Problem with camera reverse labels. Stopping here...')
                
                dataset_ok = 0;
            end
        end
        
        % add documentation
        ts.help = ['Structure designed to hold data contained in timestamp.dat from two cameras labeled cam0 and cam1. \n' ...
            'cam0/1.frames, cam0/1.times, cam0/1.N, cam0/1.buffer and cam0/1.label contain the number of frames (.frames), \n' ...
            'system clock time of each frame acquisition (.times), the total number of frames (.N) the number of frames held \n' ...
            'in buffer (.buffer), and the label of each camera (i.e. , .label = miniscope or behavior). \n' ...
            'timestamp.miniscope/behavior.label store these same labels but from the miniscope/behavior \n' ...
            'cameras point of view (i.e label = 0 or 1).'];
        
        % export the value
        newts = ts;
        
    end


% plot frames vs times
    function disp_time_vs_frames(ts,msNfr,bhNfr)
        
        % plot frames vs time
        % miniscope = red, behavior = green
        
        if isequal([ts.cam0.N ts.cam1.N],[sum(msNfr) sum(bhNfr)])
            cam_colors = ['r';'g';'r';'g'];
            color0 = [1 0 0];
            color1 = [0 1 0];
        elseif isequal([ts.cam1.N ts.cam0.N],[sum(msNfr) sum(bhNfr)])
            cam_colors = ['g';'r';'r';'g'];
            color0 = [0 1 0];
            color1 = [1 0 0];
        end
        
        % camera 0
        subplot(2,4,5)
        hold on
        
        plot(ts.cam0.frames,ts.cam0.times,'o-','Color',color0);
        
        % fit linear relationship
        M0 = [ts.cam0.frames ones(ts.cam0.N,1)];
        Minv0 = pinv(M0);
        betas0 = Minv0*ts.cam0.times;
        residuals0 = ts.cam0.times - M0*betas0;
        
        % t-values to evaluate the quality of the linear fit
        t0 = betas0/std(residuals0);
        % plot the fit
        plot(ts.cam0.frames,M0*betas0,'-','Color',color0/2);
        % and the buffer
        plot(ts.cam0.frames,ts.cam0.buffer,'-','Color',color0+[0 0 1])
        ylim([0 max(ts.cam0.times)])
        
        ylabel(['time (ms), max = ' num2str(ts.cam0.times(end)/1000) 's']);
        if strcmp(ts.cam0.label,'miniscope')
            xlabel(['Frames, miniscope (red)/ buffer (mauve) - t = ' num2str(t0(1))]);
        else
            xlabel(['Frames, behavior (green)/ buffer (cyan) - t = ' num2str(t0(1))]);
        end
        
        % camera 1
        subplot(2,4,6)
        hold on
        
        plot(ts.cam1.frames,ts.cam1.times,'o-','Color',color1);
        
        % fit linear relationship
        M1 = [ts.cam1.frames ones(ts.cam1.N,1)];
        Minv1 = pinv(M1);
        betas1 = Minv1*ts.cam1.times;
        residuals1 = ts.cam1.times - M1*betas1;
        
        % t-values to evaluate the quality of the linear fit
        t1 = betas1/std(residuals1);
        
        plot(ts.cam1.frames,M1*betas1,'-','Color',color1/2);
        % and the buffer
        plot(ts.cam1.frames,ts.cam1.buffer,'-','Color',color1+[0 0 1])
        ylim([0 max(ts.cam1.times)])
        
        ylabel(['time (ms), max = ' num2str(ts.cam1.times(end)/1000) 's']);
        if strcmp(ts.cam1.label,'miniscope')
            xlabel(['Frames, miniscope (red)/ buffer (mauve) - t = ' num2str(t1(1))]);
        else
            xlabel(['Frames, behavior (green)/ buffer (cyan) - t = ' num2str(t1(1))]);
        end
        
        
        % plot the residuals
        subplot(2,4,2)
        hold on
        plot(ts.cam0.frames,residuals0,'-o','Color',color0);
        % and the buffer
        plot(ts.cam0.frames,ts.cam0.buffer,'-','Color',color0+[0 0 1])
        
        plot(ts.cam1.frames,residuals1,'-o','Color',color1);
        % and the buffer
        plot(ts.cam1.frames,ts.cam1.buffer,'-','Color',color1+[0 0 1])
        
        xlabel('frames');
        ylabel(['residuals, miniscope (red)/ behavior (green)']);
        
        temp = union(residuals0,residuals1);
        ylim([min(temp) max(temp)])
    end


% synchronization
    function [sn] = synchronize(ts)
        
        % 1) miniscope = master, behavior = slave
        my_label = ['cam' ts.miniscope.label];
        % miniscope frames
        eval(['sn.miniscopeMaster.masterFrames = ts.' my_label '.frames;']);
        % miniscope times
        eval(['sn.miniscopeMaster.masterTimes = ts.' my_label '.times;']);
        % number of frames
        sn.miniscopeMaster.N = length(sn.miniscopeMaster.masterFrames);
        
        % synchronized behavior frames
        my_label = ['cam' ts.behavior.label];
        % time of the behavior frames
        behav_times = [];
        eval(['behav_times = ts.' my_label '.times;']);
        
        for f=1:sn.miniscopeMaster.N
            
            % time of the miniscope frame of interest
            ref_time = sn.miniscopeMaster.masterTimes(f);
            
            % look for closest behavior time
            my_diff = abs(ref_time - behav_times);
            [vmin,imin] = min(my_diff);
            
            % set it as "synchronized" frame for behavior
            if length(vmin)==1
                % put the position of the minimum as frame for behavior
                sn.miniscopeMaster.slaveFrames(f) = imin;
                % also write down the slave time for that frame
                sn.miniscopeMaster.slaveTimes(f) = behav_times(imin);
            else
                disp('Problem synchronizing behavior to miniscope: more than one minimum found');
                
                dataset_ok = 0;
            end
            
        end
        
        
        % 2) behavior =  master, miniscope =slave
        my_label = ['cam' ts.behavior.label];
        % miniscope frames
        eval(['sn.behaviorMaster.masterFrames = timestamp.' my_label '.frames;']);
        % miniscope times
        eval(['sn.behaviorMaster.masterTimes = timestamp.' my_label '.times;']);
        % number of frames
        sn.behaviorMaster.N = length(sn.behaviorMaster.masterFrames);
        
        % synchronized miniscope frames
        my_label = ['cam' ts.miniscope.label];
        % time of the miniscope frames
        mini_times = [];
        eval(['mini_times = ts.' my_label '.times;']);
        
        for f=1:sn.behaviorMaster.N
            
            % time of the behavior frame of interest
            ref_time = sn.behaviorMaster.masterTimes(f);
            
            % look for closest miniscope time
            my_diff = abs(ref_time - mini_times);
            [vmin,imin] = min(my_diff);
            
            % set it as "synchronized" frame for miniscope
            if length(vmin)==1
                % put the position of the minimum as frame for miniscope
                sn.behaviorMaster.slaveFrames(f) = imin;
                % also write down the slave time for that frame
                sn.behaviorMaster.slaveTimes(f) = mini_times(imin);
            else
                disp('Problem synchronizing miniscope to behavior: more than one minimum found');
                
                dataset_ok = 0;
            end
            
        end
        
        % information
        sn.help = ...
            ['synchronization contains 2 matrices, i.e. lookup tables, for synchronizing miniscope frames, behavior frames and cage events. \n' ...
            'Each matrix performs synchronization with a different reference time series: miniscopeMaster/behaviorMaster = synchronization \n' ...
            'using miniscope/behavior frames as reference. Each of them features .masterTimes/Frames and .slaveTimes/Frames fields \n' ...
            'with master referring to the reference, and slave to the synchronized, time series. The total number of master \n' ...
            'frames is stored in .N, and the corresponding cage times in .cage.'];
        
    end


% add cage time to synchronization matrices
    function [newsn] = compute_cage_time(sn,fschedule)
        
        % First look for evnt_name = "Output On Event" and item_name = "TTL #1"
        % in the schedule csv file, and look up its timing in the cage coordinate system.
        
        fid = fopen(fschedule);
        tline = fgetl(fid);
        TTL_time = -Inf;
        nb_lines = 0;
        
        disp(['Starting to read file ' fschedule]);
        
        while ischar(tline)
            
            nb_lines = nb_lines + 1;
            
            % read each line in full and then partition using ',' as delimiter (as
            % this is a csv file).
            tline = fgetl(fid);
            
            if not(isequal(tline,-1))
                
                % analyse the content of each line. "," is the separator
                separator_pos = find(tline==',');
                
                % modified >>>
                % first column contains the time stamp
                time = tline(1:separator_pos(1));
                time(find(time=='"')) = [];
                time = str2num(time);
                
                % fourth column contains a lot of the event stamps.
                
                % read the item name
                item_name = tline(separator_pos(3)+1:separator_pos(4)-1);
                item_name(find(item_name=='"')) = [];
                
                % read the event name
                evnt_name = tline(separator_pos(2)+1:separator_pos(3)-1);
                evnt_name(find(evnt_name=='"')) = [];
                
                % read alias name
                alias_name = tline(separator_pos(4)+1:separator_pos(5)-1);
                alias_name(find(alias_name=='"')) = [];
                
                % read arg1 value (contains the sample, target, nose poked position,
                % etc.)
                arg1_value = tline(separator_pos(8)+1:separator_pos(9)-1);
                arg1_value(find(arg1_value=='"')) = [];
                
                % look for TTL event
                if strcmp(evnt_name,'Output On Event') && strcmp(item_name,'TTL #1')
                    TTL_time = time;
                end
            end
        end
        
        fclose(fid);
        
        disp(['Done reading file ' fschedule]);
        
        if isinf(TTL_time)
            disp('No TTL event was found in the csv schedule file. Stopping here...');
            
            dataset_ok = 0;
        end
        
        % since the TTL signal triggers the behavioral and miniscope cameras, frame
        % 1 then takes place at the TTL time for the cage. So, we simply add to the
        % master clock TTL_time (in seconds).
        sn.miniscopeMaster.cageTimes = ...
            sn.miniscopeMaster.masterTimes/1000 + TTL_time - sn.miniscopeMaster.masterTimes(1)/1000;
        
        sn.behaviorMaster.cageTimes = ...
            sn.behaviorMaster.masterTimes/1000 + TTL_time - sn.behaviorMaster.masterTimes(1)/1000;
        
        % output results
        newsn = sn;
        
    end


% cleanup end of synchronization matrices
    function [newsn] = cleanup(sn)
        
        % first, miniscope = master, remove multiple identical frames in
        % slave column
        cleaned = 0;
        last_frame = sn.miniscopeMaster.slaveFrames(end);
        
        counter = sn.miniscopeMaster.N-1;
        while sn.miniscopeMaster.slaveFrames(counter)==last_frame
            % remove frame before counter
            sn.miniscopeMaster.slaveFrames(counter+1) = [];
            sn.miniscopeMaster.slaveTimes(counter+1) = [];
            sn.miniscopeMaster.masterFrames(counter+1) = [];
            sn.miniscopeMaster.masterTimes(counter+1) = [];
            sn.miniscopeMaster.cageTimes(counter+1) = [];
            sn.miniscopeMaster.N = length(sn.miniscopeMaster.masterTimes);
            
            % move on to the previous entry
            counter = counter - 1;
            cleaned = cleaned + 1;
        end
        
        if cleaned>0
            disp([num2str(cleaned) ' entries were removed from the synchronization.miniscopeMaster matrix.']);
        end
        
        % second, behavior = master, remove multiple identical frames in
        % slave column
        cleaned = 0;
        last_frame = sn.behaviorMaster.slaveFrames(end);
        
        counter = sn.behaviorMaster.N-1;
        while sn.behaviorMaster.slaveFrames(counter)==last_frame
            % remove frame before counter
            sn.behaviorMaster.slaveFrames(counter+1) = [];
            sn.behaviorMaster.slaveTimes(counter+1) = [];
            sn.behaviorMaster.masterFrames(counter+1) = [];
            sn.behaviorMaster.masterTimes(counter+1) = [];
            sn.behaviorMaster.cageTimes(counter+1) = [];
            sn.behaviorMaster.N = length(sn.behaviorMaster.masterTimes);
            
            % move on to the previous entry
            counter = counter - 1;
            cleaned = cleaned + 1;
        end
        
        if cleaned>0
            disp([num2str(cleaned) ' entries were removed from the synchronization.behaviorMaster matrix.']);
        end
        
        newsn = sn;
    end


% visual summary of both synchronization matrices
    function disp_synch_matrices(sn,ts)
        
        %  master = miniscope, slave = behavior ----------------
        subplot(2,4,[3 7])
        title(['Data Processed is in folder ' data_dir],'Interpreter','None');
        
        % master = left, slave = right
        
        % We want to normalize vertical axis (i.e. time) in [0,1]. We need
        % maximum time in the matrix
        maxT = max([max(sn.miniscopeMaster.masterTimes) max(sn.miniscopeMaster.slaveTimes)]);
        
        xlim([-1 1])
        ylim([-0.05 1.05])
        
        % designed to be wide enough to span about 10 frames at the
        % beginning and end of the duration of the experiment.
        a = 10*mean([sn.miniscopeMaster.masterTimes(end)/sn.miniscopeMaster.N ...
            sn.miniscopeMaster.slaveTimes(end)/sn.miniscopeMaster.N])/maxT;
        
        % plot interval
        b = 0.4;
        
        % Note that we display the matrix so that line number increases
        % downward.
        for i1=1:sn.miniscopeMaster.N
            % master time
            if sn.miniscopeMaster.masterTimes(i1)*[1 1]/maxT<=a
                % first lines
                line([-0.6 0.4],-b/a*sn.miniscopeMaster.masterTimes(i1)*[1 1]/maxT+1,'Color',[1 0 0]);
                text(-0.75,-b/a*sn.miniscopeMaster.masterTimes(i1)/maxT+1,num2str(sn.miniscopeMaster.masterFrames(i1)),'LineWidth',2);
            end
            if sn.miniscopeMaster.masterTimes(i1)*[1 1]/maxT>=(1-a)
                % last lines
                line([-0.6 0.4],-b/a*(sn.miniscopeMaster.masterTimes(i1)*[1 1]/maxT-1),'Color',[1 0 0]);
                text(-0.75,-b/a*(sn.miniscopeMaster.masterTimes(i1)/maxT-1),num2str(sn.miniscopeMaster.masterFrames(i1)),'LineWidth',2);
            end
            
            % slave time
            % first lines
            if sn.miniscopeMaster.slaveTimes(i1)*[1 1]/maxT<=a
                line([-0.4 0.6],-b/a*sn.miniscopeMaster.slaveTimes(i1)*[1 1]/maxT+1,'Color',[0 0 1]);
                text(0.75,-b/a*sn.miniscopeMaster.slaveTimes(i1)/maxT+1,num2str(sn.miniscopeMaster.slaveFrames(i1)),'LineWidth',2);
            end
            % last lines
            if sn.miniscopeMaster.slaveTimes(i1)*[1 1]/maxT>=(1-a)
                line([-0.4 0.6],-b/a*(sn.miniscopeMaster.slaveTimes(i1)*[1 1]/maxT-1),'Color',[0 0 1]);
                text(0.75,-b/a*(sn.miniscopeMaster.slaveTimes(i1)/maxT-1),num2str(sn.miniscopeMaster.slaveFrames(i1)),'LineWidth',2);
            end
            
            % lines joining associated frames
            if sn.miniscopeMaster.masterTimes(i1)*[1 1]/maxT<=a
                line([-0.75 0.75],[-b/a*sn.miniscopeMaster.masterTimes(i1)/maxT+1 ...
                    -b/a*sn.miniscopeMaster.slaveTimes(i1)/maxT+1],'Color',[0 1 0],'LineStyle','--');
            end
            if sn.miniscopeMaster.masterTimes(i1)*[1 1]/maxT>=(1-a)
                line([-0.75 0.75],-b/a*[sn.miniscopeMaster.masterTimes(i1)/maxT-1 ...
                    sn.miniscopeMaster.slaveTimes(i1)/maxT-1],'Color',[0 1 0],'LineStyle','--');
            end
            
        end
        
        % add cyan lines for the slave times unused in the synchronization
        if strcmp(ts.cam0.label,'miniscope')
            % miniscope is master so cam 0 is master and cam 1 is slave
            all_slave_times = ts.cam1.times;
        else
            % miniscope is master so cam 1 is master and cam 0 is slave
            all_slave_times = ts.cam0.times;
        end
        
        for i1=1:length(all_slave_times)
            if all_slave_times(i1)/maxT<=a
                if not(ismember(all_slave_times(i1),sn.miniscopeMaster.slaveTimes))
                    line([-0.4 0.6],-b/a*all_slave_times(i1)*[1 1]/maxT+1,'Color',[0 1 1])
                end
            end
            
            if all_slave_times(i1)/maxT>=1-a
                if not(ismember(all_slave_times(i1),sn.miniscopeMaster.slaveTimes))
                    line([-0.4 0.6],-b/a*(all_slave_times(i1)*[1 1]/maxT-1),'Color',[0 1 1])
                end
            end
        end
        
        set(gca,'XTick',[-0.75 0.75],'XTickLabel',{'miniscope','behavior'})
        xlabel('synchronization.miniscopeMaster')
        set(gca,'YTick',[])
        
        
        % master = behavior, slave = miniscope -----------
        subplot(2,4,[4 8])
        
        % master = right, slave = left
        
        % We want to normalize vertical axis (i.e. time) in [0,1]. We need
        % maximum time in the matrix
        maxT = max([max(sn.behaviorMaster.masterTimes) max(sn.behaviorMaster.slaveTimes)]);
        
        xlim([-1 1])
        ylim([-0.05 1.05])
        
        % designed to be wide enough to span about 10 frames at the
        % beginning and end of the duration of the experiment.
        a = 10*mean([sn.behaviorMaster.masterTimes(end)/sn.behaviorMaster.N ...
            sn.behaviorMaster.slaveTimes(end)/sn.behaviorMaster.N])/maxT;
        
        % plot interval
        b = 0.4;
        
        % Note that we display the matrix so that line number increases
        % downward.
        for i1=1:sn.behaviorMaster.N
            % slave time
            if sn.behaviorMaster.slaveTimes(i1)*[1 1]/maxT<=a
                % first lines
                line([-0.6 0.4],-b/a*sn.behaviorMaster.slaveTimes(i1)*[1 1]/maxT+1,'Color',[1 0 0]);
                text(-0.75,-b/a*sn.behaviorMaster.slaveTimes(i1)/maxT+1,num2str(sn.behaviorMaster.slaveFrames(i1)),'LineWidth',2);
            end
            if sn.behaviorMaster.slaveTimes(i1)*[1 1]/maxT>=(1-a)
                % last lines
                line([-0.6 0.4],-b/a*(sn.behaviorMaster.slaveTimes(i1)*[1 1]/maxT-1),'Color',[1 0 0]);
                text(-0.75,-b/a*(sn.behaviorMaster.slaveTimes(i1)/maxT-1),num2str(sn.behaviorMaster.slaveFrames(i1)),'LineWidth',2);
            end
            
            % master time
            % first lines
            if sn.behaviorMaster.masterTimes(i1)*[1 1]/maxT<=a
                line([-0.4 0.6],-b/a*sn.behaviorMaster.masterTimes(i1)*[1 1]/maxT+1,'Color',[0 0 1]);
                text(0.75,-b/a*sn.behaviorMaster.masterTimes(i1)/maxT+1,num2str(sn.behaviorMaster.masterFrames(i1)),'LineWidth',2);
            end
            % last lines
            if sn.behaviorMaster.masterTimes(i1)*[1 1]/maxT>=(1-a)
                line([-0.4 0.6],-b/a*(sn.behaviorMaster.masterTimes(i1)*[1 1]/maxT-1),'Color',[0 0 1]);
                text(0.75,-b/a*(sn.behaviorMaster.masterTimes(i1)/maxT-1),num2str(sn.behaviorMaster.masterFrames(i1)),'LineWidth',2);
            end
            
            % lines joining associated frames
            if sn.behaviorMaster.masterTimes(i1)*[1 1]/maxT<=a
                line([-0.75 0.75],[-b/a*sn.behaviorMaster.slaveTimes(i1)/maxT+1 ...
                    -b/a*sn.behaviorMaster.masterTimes(i1)/maxT+1],'Color',[0 1 0],'LineStyle','--');
            end
            if sn.behaviorMaster.masterTimes(i1)*[1 1]/maxT>=(1-a)
                line([-0.75 0.75],-b/a*[sn.behaviorMaster.slaveTimes(i1)/maxT-1 ...
                    sn.behaviorMaster.masterTimes(i1)/maxT-1],'Color',[0 1 0],'LineStyle','--');
            end
            
        end
        
        % add cyan lines for the slave times unused in the synchronization
        if strcmp(ts.cam0.label,'behavior')
            % behavior is master so cam 0 is master and cam 1 is slave
            all_slave_times = ts.cam1.times;
        else
            % behavior is master so cam 1 is master and cam 0 is slave
            all_slave_times = ts.cam0.times;
        end
        for i1=1:length(all_slave_times)
            % first few lines
            if all_slave_times(i1)/maxT<=a
                if not(ismember(all_slave_times(i1),sn.behaviorMaster.slaveTimes))
                    line([-0.4 0.6],-b/a*all_slave_times(i1)*[1 1]/maxT+1,'Color',[0 1 1])
                end
            end
            % last few lines
            if all_slave_times(i1)/maxT>=1-a
                if not(ismember(all_slave_times(i1),sn.behaviorMaster.slaveTimes))
                    line([-0.4 0.6],-b/a*(all_slave_times(i1)*[1 1]/maxT-1),'Color',[0 1 1])
                end
            end
        end
        
        set(gca,'XTick',[-0.75 0.75],'XTickLabel',{'miniscope','behavior'})
        xlabel('synchronization.behaviorMaster')
        set(gca,'YTick',[])
        
    end

% add "events" and save synchronization to a .mat file
    function [newsn] = add_events(sn,fschedule)
        
        % read the cage output file in .csv format and extract its
        % content in arrays
        eventTime = zeros(1,nb_lines-1);   % Touchscreen event timestamps
        eventname = strings(1,nb_lines-1);
        item_name = strings(1,nb_lines-1);  % Touchscreen event name
        aliasname = strings(1,nb_lines-1);
        groupID = zeros(1,nb_lines-1);  %Touchscreen event group ID
        arg = strings(1,nb_lines-1);
        
        fid = fopen(fschedule,'r');
        metanames = textscan(fid,'%s%s%f','delimiter',',');
        fclose(fid);
        
        Tschedule = readcell(fschedule);
        
        eventTime = cell2mat(Tschedule(length(metanames)+1:end, 1)); %time is in: ___
        eventname = string(Tschedule(length(metanames)+1:end, 3));        
        item_name = string(Tschedule(length(metanames)+1:end, 4));
        aliasname = string(Tschedule(length(metanames)+1:end, 5)); %What is alias name?
        groupID = cell2mat(Tschedule(length(metanames)+1:end, 6));
        arg = cell2mat(Tschedule(length(metanames)+1:end,9));                                               
        
        % behavior frames and times
        frameMap = sn.miniscopeMaster.slaveFrames;
        timeMap = sn.miniscopeMaster.slaveTimes/1000;
        frameMap = frameMap';
        timeMap = timeMap';
        
        % used below to synchronize the cage events to the miniscope
        a = [frameMap sn.miniscopeMaster.masterTimes sn.miniscopeMaster.slaveTimes'];
        
        loc = strfind(item_name,'TTL');                                             % find the TTL initialization
        startTime = NaN;
        startInd = NaN;
        endTime = NaN;
        endInd = NaN;
        firstTTL = 0;
        
        % find miniscope start point in behaviour
        for t = 1 : length(eventTime)
            if ~isempty(loc{t}) && firstTTL ==0
                startTime = eventTime(t);
                startInd = t;
                firstTTL = 1;
            elseif ~isempty(loc{t}) && firstTTL >0
                endTime = eventTime(t);
                endInd = t;
            end
        end
        
        eventInd = [];
        timeMap(:,1)  = timeMap(:,1) + startTime;                                   % Add the miniscope start time to miniscope time
        
        % behaviour event time to miniscope time
        for t = 1 : length(eventTime)
            if eventTime(t) <= max(timeMap)
                [minDiff, ind ] = min(abs(eventTime(t) - timeMap));
                eventInd(t) = ind;
            end
        end
        eventInd = eventInd';
        
        my_events = strings(length(frameMap),6);
        my_events(:,1) = string(frameMap);
        my_events(eventInd(startInd),2) = timeMap(eventInd(startInd));
        my_events(eventInd(startInd),3) = item_name(startInd);
        ind = 1;
        
        
        % Then, we fill the columns for the cage events
        for i = startInd : length(eventInd)
            
            if abs(ind - eventInd(i))>0
                
                group = find(eventInd == eventInd(i));
                [c,d] = strtok(item_name(group(find(contains(item_name(group),'Correction_Trial'),1))),'Correction_Trials');
                [a1,a2]=strtok(arg(group(find(contains(item_name(group),'Correction_Trial'),1))),'Value');
                
                if ~isempty(find(contains(item_name(group),'BIRBeam'),1)) && ~isempty(find(contains(item_name(group),'Correction_Trial'),1))&& c=='' && d=='' && a1 == 1 && ~isempty(find(contains(eventname(group),'Display Image'),1)) && ~isempty(find(contains(aliasname(group),'Training'),1))
                    my_events(eventInd(i),2) = timeMap(eventInd(i));
                    my_events(eventInd(i),3) = 'Display Image & BIRBeam#1 & correction Trial Start';
                    [a,b]=strtok(arg(group(find(contains(eventname(group),'Display Image'),1))),['Image' 'Position']);
                    my_events(eventInd(i),4) = strtok(a);
                elseif ~isempty(find(contains(eventname(group),'Display Image'),1)) && ~isempty(find(contains(aliasname(group),'Training'),1)) && ~isempty(find(contains(item_name(group),'BIRBeam'),1))
                    my_events(eventInd(i),2) = timeMap(eventInd(i));
                    my_events(eventInd(i),3) = 'BIRBeam & Display Image';
                    [a,b]=strtok(arg(group(find(contains(eventname(group),'Display Image'),1))),['Image' 'Position']);
                    my_events(eventInd(i),4) = strtok(a);
                elseif ~isempty(find(contains(eventname(group),'Display Image'),1)) && ~isempty(find(contains(aliasname(group),'Training'),1))
                    my_events(eventInd(i),2) = timeMap(eventInd(i));
                    my_events(eventInd(i),3) = 'Display Image';
                    [a,b]=strtok(arg(group(find(contains(eventname(group),'Display Image'),1))),['Image' 'Position']);
                    my_events(eventInd(i),4) = strtok(a);
                elseif ~isempty(find(contains(item_name(group),'BIRBeam'),1)) && ~isempty(find(contains(item_name(group),'Correction_Trial'),1)) && a1 == 1 && c=='' && d==''
                    my_events(eventInd(i),2) = timeMap(eventInd(i));
                    my_events(eventInd(i),3) = 'BIRBeam#1 & correction Trial Start';
                elseif ~isempty(find(contains(item_name(group),'BIRBeam'),1))
                    my_events(eventInd(i),2) = timeMap(eventInd(i));
                    my_events(eventInd(i),3) = 'BIRBeam#1';
                elseif ~isempty(find(contains(eventname(group),'Touch Down Event'),1)) && ~isempty(find(contains(item_name(group),'Incorrect'),1)) && ~isempty(find(contains(item_name(group),'Start Delay'),1))
                    my_events(eventInd(i),2) = timeMap(eventInd(i));
                    my_events(eventInd(i),3) = 'Nose-Poke Incorrect';
                    [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);
                    my_events(eventInd(i),4) = strtok(a);
                    my_events(eventInd(i),6) = 1;
                elseif ~isempty(find(contains(eventname(group),'Touch Down Event'),1)) && ~isempty(find(contains(item_name(group),'Correct'),1)) && ~isempty(find(contains(item_name(group),'Start Delay'),1))
                    my_events(eventInd(i),2) = timeMap(eventInd(i));
                    my_events(eventInd(i),3) = 'Nose-Poke Correct';
                    [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);
                    my_events(eventInd(i),4) = strtok(a);
                    my_events(eventInd(i),6) = 1;
                elseif ~isempty(find(contains(eventname(group),'Touch Down Event'),1)) && ~isempty(find(contains(item_name(group),'Incorrect'),1))
                    my_events(eventInd(i),2) = timeMap(eventInd(i));
                    my_events(eventInd(i),3) = 'Nose-Poke Incorrect';
                    [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);
                    my_events(eventInd(i),4) = strtok(a);
                elseif ~isempty(find(contains(eventname(group),'Touch Down Event'),1)) && ~isempty(find(contains(item_name(group),'Correct'),1))
                    my_events(eventInd(i),2) = timeMap(eventInd(i));
                    my_events(eventInd(i),3) = 'Nose-Poke Correct';
                    [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);
                    my_events(eventInd(i),4) = strtok(a);
                elseif ~isempty(find(contains(eventname(group),'Touch Down Event'),1))
                    if max(group) < length(eventInd)
                        group2 = find(eventInd == eventInd(max(group)+1));
                    else
                        group2 = [];
                    end
                    if ~isempty(group2) && ~isempty(find(strcmp(item_name(group2),'Incorrect'),1)) && ~isempty(find(contains(item_name(group),'Start Delay'),1))
                        my_events(eventInd(i),2) = timeMap(eventInd(i));
                        my_events(eventInd(i),3) = 'Nose-Poke Incorrect';
                        [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);
                        my_events(eventInd(i),4) = strtok(a);
                        my_events(eventInd(i),6) = 1;
                    elseif ~isempty(group2) && ~isempty(find(strcmp(item_name(group2),'Correct'),1))&& ~isempty(find(contains(item_name(group),'Start Delay'),1))
                        my_events(eventInd(i),2) = timeMap(eventInd(i));
                        my_events(eventInd(i),3) = 'Nose-Poke Correct';
                        [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);
                        my_events(eventInd(i),4) = strtok(a);
                        my_events(eventInd(i),6) = 1;
                    elseif ~isempty(group2) && ~isempty(find(strcmp(item_name(group2),'Incorrect'),1)) && isempty(find(contains(item_name(group2),'Next trial'),1))
                        my_events(eventInd(i),2) = timeMap(eventInd(i));
                        my_events(eventInd(i),3) = 'Nose-Poke Incorrect';
                        [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);
                        my_events(eventInd(i),4) = strtok(a);
                    elseif ~isempty(group2) && ~isempty(find(strcmp(item_name(group2),'Correct'),1)) && isempty(find(contains(item_name(group2),'Next trial'),1))
                        my_events(eventInd(i),2) = timeMap(eventInd(i));
                        my_events(eventInd(i),3) = 'Nose-Poke Correct';
                        [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);
                        my_events(eventInd(i),4) = strtok(a);
                    elseif ~isempty(find(contains(item_name(group),'Start Delay'),1))
                        my_events(eventInd(i),2) = timeMap(eventInd(i));
                        my_events(eventInd(i),3) = 'Nose-Poke';
                        [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);
                        my_events(eventInd(i),4) = strtok(a);
                        my_events(eventInd(i),6) = 1;
                    else
                        my_events(eventInd(i),2) = timeMap(eventInd(i));
                        my_events(eventInd(i),3) = 'Nose-Poke';
                        [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Down Event')+i-1,1))),['Image' 'Position']);
                        my_events(eventInd(i),4) = strtok(a);
                    end
                elseif ~isempty(find(contains(eventname(group),'Touch Up'),1))
                    my_events(eventInd(i),2) = timeMap(eventInd(i));
                    my_events(eventInd(i),3) = eventname(i);
                    [a,b]=strtok(arg(group(find(contains(eventname(group),'Touch Up Event')+i-1,1))),['Image' 'Position']);
                    my_events(eventInd(i),4) = strtok(a);
                else
                    my_events(eventInd(i),2) = timeMap(eventInd(i));
                    my_events(eventInd(i),3) = item_name(i);
                end
                if ~isempty(find(contains(item_name(group),'Correction_Trial'),1))&& a1 == 1 && c=='' && d==''
                    my_events(eventInd(i),2) = timeMap(eventInd(i));
                    my_events(eventInd(i),3) = 'correction Trial';
                elseif ~isempty(find(contains(item_name(group),'Start Delay'),1))
                    my_events(eventInd(i),2) = timeMap(eventInd(i));
                    my_events(eventInd(i),6) = 1;
                elseif ~isempty(find(contains(item_name(group),'Delay End'),1))
                    my_events(eventInd(i),2) = timeMap(eventInd(i));
                    my_events(eventInd(i),6) = 2;
                end
                if ~isempty(find(contains(item_name(group),'Next trial'),1)) &&   strcmp(item_name(group(find(contains(item_name(group),'Next trial'),1))), 'Next trial')
                    my_events(eventInd(i),5) = '1';
                end
            end
        end
        
        % finally, add events
        sn.miniscopeMaster.events = my_events;
        
        % output the result
        newsn = sn;
    end

% save figures
    function save_figures()
        
        figure(h_synchr);
        img = getframe(gcf);
        imwrite(img.cdata,[data_dir '\summary_synchronization.png']);
        savefig([data_dir '\summary_synchronization.fig']);
        
        figure(h_videos);
        img = getframe(gcf);
        imwrite(img.cdata,[data_dir '\summary_videos.png']);
        savefig([data_dir '\summary_videos.fig']);
        
        figure(h_timestamp);
        img = getframe(gcf);
        imwrite(img.cdata,[data_dir '\summary_timestamp.png']);
        savefig([data_dir '\summary_timestamp.fig']);
        
    end

end