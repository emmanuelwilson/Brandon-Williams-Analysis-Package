
%**************************************************************************
% Only part you need to change. The trajectory and bad data marker are
% saved where the data is, and also a couple of summary images of the
% data/tracking quality.
% *************************************************************************

%paths{1} = 'D:\Andres\100203\TUNL\Stage 1\S3\6 SEC\100203_10feb19_H10_M44_S17\';

function [] = extract_clean_func_crawl(p,visual_check)
% visual check of the extraction
% visual_check = 0;
paths = genpath(p);
paths = strsplit(paths,';')';

% scale factor to convert pixels into mm in the cage
fact = 238/450;

for directory = 1:length(paths)
    if ~isempty(paths{directory})
        disp(' ');
        disp('--------------------------------------------------------------');
        my_path = paths{directory};
        
        disp([num2str(directory) '/: ' my_path]);
        disp(' ');
        
        d = dir(paths{directory});
        fnames = {d.name};
        if isempty(find(contains(fnames,'behavCam'),1))
            disp('No Behaviour Videos Found ');
        else                        
            all_max = [];
            mean_int = [];
            
            my_lines = [];
            my_columns = [];
            
            % count the number of videos in the folder
            temp=dir(my_path);
            vid_nums = 0;
            for i1=1:length(temp)
                if length(temp(i1).name)>10
                    temp1 = temp(i1).name;
                    if strcmp(temp1(1:8),'behavCam')==1 && strcmp(temp1(end-3:end),'.avi')==1
                        vid_nums = vid_nums+1;
                    end
                end
            end
            
            Nvid = vid_nums;
            total_frames = 0;
            
            % detect the position of the image maximum for each frame
            for v=1:Nvid
                
                my_file = [my_path '/behavCam' num2str(v) '.avi'];
                
                disp([num2str(v) '/' num2str(Nvid) ' - ' my_file]);
                
                % load the video and convert to double format.
                obj = VideoReader(my_file);
                one_movie = double(obj.read());
                
                % look for the position of the maximum
                bw = squeeze(mean(one_movie,3));
                Nfr = size(bw,3);
                
                total_frames = total_frames + Nfr;
                
                for f=1:Nfr
                    temp = bw(:,:,f);
                    mean_int = [mean_int mean(temp(:))];
                    mmax = max(temp(:));
                    [l,c] = ind2sub([size(bw,1),size(bw,2)],find(mmax==bw(:,:,f)));
                    
                    my_lines = [my_lines mean(l)];
                    my_columns = [my_columns mean(c)];
                    
                    all_max = [all_max mmax];
                end
                
            end
            
            % number of lines and columns of the video
            Nl = size(one_movie,1);
            Nc = size(one_movie,2);
            
            % threshold the variables computed above so we can separate cases where
            % the led is visible from when it is missing.
            % For this, we use the histogram of the mean intensity: it has 4 peaks.
            % One for the dark, one for sample (1 pad on), one for choice period (2
            % pads on) and one for the houselight.
            [u,v] = hist(mean_int,100);
            % u = frequency, v = mean intensity value
            
            % sorting this
            temp = [v' u'];
            [x,y] = sort(temp(:,2),'descend');
            stemp = temp(y,:);
            
            % locating when house light is on
            houselight = +(mean_int>4*stemp(1,1));
            
            % and when the LED is visible
            visible = +(all_max>200);
            
            % divide in blocks of 0's bordered by 1's
            [labeled_blocks,numRegions] = bwlabel(1-visible,4);
            
            % store the position of the led
            trajectory = [];
            trajectory.x = zeros(1,total_frames);
            trajectory.y = zeros(1,total_frames);
            
            % put the actual position when the led is visible and NaN otherwise
            for t=1:total_frames
                if labeled_blocks(t)==0
                    % x has to be largest for pad 5 (which corresponds to column 1)
                    % used to be trajectory.x(t) = my_columns(t).
                    trajectory.x(t) = Nc - my_columns(t) + 1;
                    trajectory.y(t) = my_lines(t);
                else
                    trajectory.x(t) = NaN;
                    trajectory.y(t) = NaN;
                end
            end
            
            % linear interpolation
            % interpolated over NaNs that are in the trajectories
            % .x
            nan_traj = isnan(trajectory.x);
            traj_t = 1:numel(trajectory.x);
            interp_traj.x = trajectory.x;
            interp_traj.x(nan_traj) = round(interp1(traj_t(~nan_traj), interp_traj.x(~nan_traj), traj_t(nan_traj)));
            
            % .y
            nan_traj = isnan(trajectory.y);
            traj_t = 1:numel(trajectory.y);
            interp_traj.y = trajectory.y;
            interp_traj.y(nan_traj) = round(interp1(traj_t(~nan_traj), interp_traj.y(~nan_traj), traj_t(nan_traj)));
            
            % fix NaNs at both ends of all_traj
            % .x
            nanlistx = find(isnan(interp_traj.x));
            if not(isempty(nanlistx))
                % check if led invisible at the start
                if nanlistx(1)==1
                    temp = find(not(isnan(interp_traj.x)),1,'first');
                    interp_traj.x(1:temp-1) = interp_traj.x(temp);
                end
                
                % check if led invisible at the start
                if nanlistx(end)==total_frames
                    temp = find(not(isnan(interp_traj.x)),1,'last');
                    interp_traj.x(temp+1:end) = interp_traj.x(temp);
                end
            end
            
            % .y
            nanlisty = find(isnan(interp_traj.y));
            if not(isempty(nanlisty))
                % check if led invisible at the start
                if nanlisty(1)==1
                    temp = find(not(isnan(interp_traj.y)),1,'first');
                    interp_traj.y(1:temp-1) = interp_traj.y(temp);
                end
                
                % check if led invisible at the start
                if nanlisty(end)==total_frames
                    temp = find(not(isnan(interp_traj.y)),1,'last');
                    interp_traj.y(temp+1:end) = interp_traj.y(temp);
                end
            end
            
            % next, do a low pass filtering using filtfilt and butterworth filter
            % sample freq is about 30 fps.
            Fs = 30;
            
            % maximum of 5 outlier frames (= half wave of period equal to 10 frames then)
            limit_freq = 10/Fs;
            
            d1 = designfilt('lowpassiir','FilterOrder',12, ...
                'HalfPowerFrequency',limit_freq,'DesignMethod','butter');
            filt_traj.x = filtfilt(d1,interp_traj.x);
            filt_traj.y = filtfilt(d1,interp_traj.y);
            
            % hold on
            % plot(interp_traj.x,'-')
            % plot(y,'-')
            
            
            % finally, we want to remove from the analysis any period where the led has
            % disappeared for more than 5 consecutive frames, and distance covered is
            % larger than 10mm (i.e. 20 pixels).
            
            % compute the duration and distance spanned by each block in labeled_blocks
            % compute duration of the block where the led disappears
            durations = zeros(1,numRegions);
            for i1=1:numRegions
                durations(i1) = length(find(labeled_blocks==i1));
            end
            
            % compute the distance between the beginning and end of the block
            distances = zeros(1,numRegions);
            
            for s=1:numRegions
                temp = find(labeled_blocks==s);
                
                % led position immediately before and after the block
                initial = temp(1)-1;
                if initial==0
                    initial = 1;
                end
                final = temp(end)+1;
                if final>total_frames
                    final = total_frames;
                end
                p1 = round([my_lines(initial) my_columns(initial)]);
                p2 = round([my_lines(final) my_columns(final)]);
                my_dist = p1-p2;
                my_dist = sqrt(sum(my_dist.^2));
                my_dist = round(my_dist);
                
                distances(s) = my_dist;
            end
            
            % finally make a mask for the data we keep: durations<5 and
            % distances<10cm/fact
            bad_data = zeros(1,total_frames);
            for s=1:numRegions
                my_block = find(labeled_blocks==s);
                if durations(s)>5 || distances(s)>(100/fact)
                    bad_data(my_block) = 1;
                end
            end
            
            % also reject position when house light is on as we can't locate the
            % led when it is on
            for f=1:total_frames
                if houselight(f)==1
                    bad_data(f) = 1;
                end
            end
            
            % save the results
            save([my_path '/filt_traj.mat'],'filt_traj');
            save([my_path '/bad_data.mat'],'bad_data');
            
            % display stats on where and how much the bad trajectory is
            figure('Position',[300 500 1000 500])
            subplot(1,2,1)
            hold on
            for f=1:total_frames
                if bad_data(f)==1
                    if houselight(f)==1
                        plot(filt_traj.x(f),filt_traj.y(f),'.y')
                    else
                        plot(filt_traj.x(f),filt_traj.y(f),'.r')
                    end
                else
                    plot(filt_traj.x(f),filt_traj.y(f),'.g')
                end
            end
            
            subplot(1,2,2)
            hold on
            
            nb_houselight = length(find(houselight==1));
            nb_bad = length(find(bad_data==1)) - nb_houselight;
            nb_good = length(find(bad_data==0));
            
            temp = [nb_good nb_houselight nb_bad]/total_frames;
            bar(1,temp(1),'g')
            bar(2,temp(2),'y');
            bar(3,temp(3),'r')
            set(gca,'XTick',[1 2 3],'XTickLabel',{'included','house light','rejected'})
            title(my_path)
            drawnow
            img = getframe(gcf);
            imwrite(img.cdata,[my_path '/summary_bad_data.png']);
            close all
            
            % also plot the mean intensity, the houselight extracted from it, and
            % the bad data mask
            figure
            hold on
            plot(mean_int)
            hold on
            plot(houselight*12)
            plot(bad_data*8)
            legend('<intensity>','house light','bad data');
            title(my_path)
            drawnow
            img = getframe(gcf);
            imwrite(img.cdata,[my_path '/summary_masking.png']);
            close all
            
            % It would be nice to have some kind of means to verify a little what
            % was produced.
            if visual_check==1
                
                figure
                frame = 0;
                for v=1:Nvid
                    
                    my_file = [my_path '/behavCam' num2str(v) '.avi'];
                    
                    disp([num2str(v) '/' num2str(Nvid) ' - ' my_file]);
                    
                    % load the video and convert to double format.
                    obj = VideoReader(my_file);
                    one_movie = obj.read();
                    
                    % display each frame with a cross indicating the believe led
                    % position
                    for t=1:size(one_movie,4)
                        frame = frame + 1;
                        temp = squeeze(single(one_movie(:,:,:,t))/255);
                        % converting (x,y) into (lines,columns)
                        L = round(filt_traj.y(frame));
                        C = Nc - round(filt_traj.x(frame)) + 1;
                        
                        if bad_data(frame)==0
                            temp(:,C,1) = 0;
                            temp(:,C,2) = 1;
                            temp(:,C,3) = 0;
                            temp(L,:,1) = 0;
                            temp(L,:,2) = 1;
                            temp(L,:,3) = 0;
                            image(temp);
                        else
                            temp(:,C,1) = 1;
                            temp(:,C,2) = 0;
                            temp(:,C,3) = 0;
                            temp(L,:,1) = 1;
                            temp(L,:,2) = 0;
                            temp(L,:,3) = 0;
                            image(temp);
                        end
                        xlabel(['Video ' num2str(v) ', frame ' num2str(t)]);
                        
                        drawnow
                    end
                    
                end
                
            end
            
        end
    end
end
end
