function ms = msNormCorre_05GRIN(ms,isnonrigid,doWindowed,varargin);
% Performs fast, rigid registration (option for non-rigid also available).
% Relies on NormCorre (Paninski lab). Rigid registration works fine for
% large lens (1-2mm) GRIN lenses, while non-rigid might work better for
% smaller lenses. Ideally you want to compare both on a small sample before
% choosing one method or the other.
% Original script by Eftychios Pnevmatikakis, edited by Guillaume Etter

warning off all

%% Auto-detect operating system
if ispc
    separator = '\'; % For pc operating systems
else
    separator = '/'; % For unix (mac, linux) operating systems
end

%% Filtering parameters
gSig = 7/ms.ds; % 7
gSiz = 17/ms.ds; % 17
psf = fspecial('gaussian', round(2*gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;
bound = round(ms.height/(2*ms.ds));

template = [];

if isempty(varargin)
    writerObj = VideoWriter([ms.dirName separator 'msvideo.avi'],'Grayscale AVI');
else
    checkP(varargin{1})
    writerObj = VideoWriter(varargin{1},'Grayscale AVI');
end
open(writerObj);

ms.shifts = [];
ms.meanFrame = [];
limit = [];
mask = [];
cleanDots = [];
for video_i = 1:ms.numFiles;
    if ~isempty(varargin) && video_i > 5
        continue        
    end
    
    name = [ms.vidObj{1, video_i}.Path separator ms.vidObj{1, video_i}.Name];
    fprintf(['\t\tRegistration on:  (' num2str(video_i) ' of ' num2str(ms.numFiles) ')  ' ms.dirName]);
    
    % read data and convert to single
    Yf = read_file(name);
    Yf = single(Yf);
    Yf = downsample_data(Yf,'space',1,ms.ds,1); % ms.ds instead of 2
    
    if doWindowed
        if isempty(limit)
            figure(1)
            imagesc(Yf(:,:,1));
            colormap gray
            limit = getrect();
        end
        Yf = Yf(nanmax(ceil(limit(2)),1):nanmin(ceil(limit(2)+limit(4)),length(Yf(:,1,1))),...
            nanmax(ceil(limit(1)),1):nanmin(ceil(limit(1)+limit(3)),length(Yf(1,:,1))),:);
    end
    bound = round(nanmin(size(Yf(:,1,1)))/2);
    Y = imfilter(Yf,psf,'symmetric');
    [d1,d2,T] = size(Y);
    
    if isempty(mask)
        figure(1)
        imagesc(Y(:,:,1));
        colormap parula
        colorbar
        mask = [50, 150, 500, 400];
    end

    
%     if doWindowed
%         if isempty(cleanDots)
%             for i = 1:2
%                 figure(1)
%                 imagesc(Y(:,:,500));
%                 colormap parula
%                 cleanDots(:,i) = getrect();
%             end
%         end
%         for i = 1:length(cleanDots(1,:))
%             Y(floor(cleanDots(2,i)):ceil(cleanDots(2,i)+cleanDots(4,i)),...
%                 floor(cleanDots(1,i)):ceil(cleanDots(1,i)+cleanDots(3,i)),:) = 0;% ...
% %                 repmat(Y(ceil(cleanDots(2,i)+cleanDots(4,i))+1,...
% %                 ceil(cleanDots(1,i)+cleanDots(3,i))+1,:),...
% %                 [length(floor(cleanDots(2,i)):ceil(cleanDots(2,i)+cleanDots(4,i))) ...
% %                 length(floor(cleanDots(1,i)):ceil(cleanDots(1,i)+cleanDots(3,i)))]);
%         end
%     end
      
    % Setting registration parameters (rigid vs non-rigid)
    if isnonrigid
%         disp('Non-rigid motion correction...');
        options = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',50, ...
            'grid_size',[128,128]*2,'mot_uf',4,'correct_bidir',false, ...
            'overlap_pre',32,'overlap_post',32,'max_shift',20,'iter',1); % max shift 20 grid size [128 128]*2 mot_uf 4
    else
%         disp('Rigid motion correction...');
    options = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',200,'max_shift',20,'iter',1,'correct_bidir',false);
    end
    
    %% register using the high pass filtered data and apply shifts to original data
    if isempty(template); 
        [M1,shifts1,template] = normcorre(Y(nanmax(ceil(mask(2)),1):nanmin(ceil(mask(2)+mask(4)),length(Y(:,1,1))),...
            nanmax(ceil(mask(1)),1):nanmin(ceil(mask(1)+mask(3)),length(Y(1,:,1))),:),options); % register filtered data
        % exclude boundaries due to high pass filtering effects
    else
        [M1,shifts1,template] = normcorre(Y(nanmax(ceil(mask(2)),1):nanmin(ceil(mask(2)+mask(4)),length(Y(:,1,1))),...
            nanmax(ceil(mask(1)),1):nanmin(ceil(mask(1)+mask(3)),length(Y(1,:,1))),:),options,template); % register filtered data
    end
    
    Mr = apply_shifts(Yf,shifts1,options,bound/2,bound/2); % apply shifts to full dataset
    % apply shifts on the whole movie    
    
    
    if size(Mr,3) == 3
        Mr = permute(Mr,[1 2 4 3]);
    end
    
    writeVideo(writerObj,uint8(Mr));
    
    if size(Mr,4) == 3
        Mr = permute(Mr,[1 2 4 3]);
    end
    
    %% compute metrics
    
    if video_i == ms.numFiles && length(size(Yf))<3
        ms.shifts{video_i} = shifts1;
        continue
    end
    
    [cYf,mYf,vYf] = motion_metrics(Yf,options.max_shift); 
    [cM1f,mM1f,vM1f] = motion_metrics(Mr,options.max_shift);
    
    if video_i == 1;
        ms.meanFrame = mM1f;
    else;
        ms.meanFrame = (ms.meanFrame + mM1f)./2;
    end
    corr_gain = cYf./cM1f*100;
    
    ms.shifts{video_i} = shifts1;
       
end

close(writerObj);

end