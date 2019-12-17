    
function stablizeCalcium(folder)
    % demo file for applying the NoRMCorre motion correction algorithm on 
    % 1-photon widefield imaging data
    % Example files can be obtained through the miniscope project page
    % www.miniscope.org

%     gcp;
    %% read data and convert to double
    
    folder = [folder '/msCam*.avi'];
    %addpath(genpath('../../NoRMCorre'));


    %%% Sort Clips in chronological order
    name = dir(folder);
    clipNum = nan(1,length(name));
    for i = 1:length(name)
        clipNum(i) = str2num(name(i).name(6:end-4));
    end
    [a b] = sort(clipNum);
    name = name(b);
    clipNum = clipNum(b);

%     
%     alreadyDone = true;
%     for i = 1:length(name)
%         if ~(exist([folder(1:end-10) 'Stablized_MAT/Normcorred_' num2str(clipNum(i)) '.mat'])==2);
%             alreadyDone = false;
%         end
%     end
%     
%     if alreadyDone
%         fprintf('\n\t********  ALREADY STABLIZED  *********\n');
%         return
%     end
    
    template1 = [];
    for iter = 1
        for i = 1:length(name)
            Yf = read_file([folder(1:end-10) name(i).name]);
            Yf = single(Yf);
            [d1,d2,T] = size(Yf);

            %% perform some sort of deblurring/high pass filtering

            if (0)    
                hLarge = fspecial('average', 40);
                hSmall = fspecial('average', 2); 
                for t = 1:T
                    Y(:,:,t) = filter2(hSmall,Yf(:,:,t)) - filter2(hLarge, Yf(:,:,t));
                end
                %Ypc = Yf - Y;
                bound = size(hLarge,1);
            else
                gSig = 7; 
                gSiz = 17; 
                psf = fspecial('gaussian', round(2*gSiz), gSig);
                ind_nonzero = (psf(:)>=max(psf(:,1)));
                psf = psf-mean(psf(ind_nonzero));
                psf(~ind_nonzero) = 0;   % only use pixels within the center disk
                %Y = imfilter(Yf,psf,'same');
                %bound = 2*ceil(gSiz/2);
                Y = imfilter(Yf,psf,'symmetric');
                bound = 0;
            end
            %% first try out rigid motion correction
                % exclude boundaries due to high pass filtering effects
            options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',200,'max_shift',30,'iter',1,'correct_bidir',false);


            %% register using the high pass filtered data and apply shifts to original data
            tic; [M1,shifts1,template1] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r,template1); toc % register filtered data
%                 exclude boundaries due to high pass filtering effects

    %         options_r.output_type = 'mat';
    %         options_r.mem_filename = [folder(1:end-10) name(i).name(1:end-4) '_stablized.mat'];

            tic; Mr = apply_shifts(Yf,shifts1,options_r,bound/2,bound/2); toc % apply shifts to full dataset            
%                 apply shifts on the whole movie
            %% compute metrics 
%             [cY,mY,vY] = motion_metrics(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r.max_shift);
%             [cYf,mYf,vYf] = motion_metrics(Yf,options_r.max_shift);
% 
%             [cM1,mM1,vM1] = motion_metrics(M1,options_r.max_shift);
%             [cM1f,mM1f,vM1f] = motion_metrics(Mr,options_r.max_shift);

            %% plot rigid shifts and metrics
%             shifts_r = squeeze(cat(3,shifts1(:).shifts));
    %         figure;
    %             subplot(311); plot(shifts_r);
    %                 title('Rigid shifts','fontsize',14,'fontweight','bold');
    %                 legend('y-shifts','x-shifts');
    %             subplot(312); plot(1:T,cY,1:T,cM1);
    %                 title('Correlation coefficients on filtered movie','fontsize',14,'fontweight','bold');
    %                 legend('raw','rigid');
    %             subplot(313); plot(1:T,cYf,1:T,cM1f);
    %                 title('Correlation coefficients on full movie','fontsize',14,'fontweight','bold');
    %                 legend('raw','rigid');

    
            %%%% NONRIGID
    
%             options_nr = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',50, ...
%                 'grid_size',[128,128]*2,'mot_uf',4,'correct_bidir',false, ...
%                 'overlap_pre',32,'overlap_post',32,'max_shift',30);
% 
%             tic; [M1,shifts2,template1] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_nr,template1); toc % register filtered data
%             tic; Mr = apply_shifts(Yf,shifts2,options_nr,bound/2,bound/2); toc % apply the shifts to the removed percentile

            %% compute metrics

%             [cM2,mM2,vM2] = motion_metrics(M2,options_nr.max_shift);
%             [cM2f,mM2f,vM2f] = motion_metrics(Mpr,options_nr.max_shift);
%     
            %% display downsampled data

            tsub = 5;

            Y_ds = downsample_data(Y,'time',tsub);
            Yf_ds = downsample_data(Yf,'time',tsub);
            M1_ds = downsample_data(M1,'time',tsub);
            M1f_ds = downsample_data(Mr,'time',tsub);
        %     M2_ds = downsample_data(M2,'time',tsub);
        %     M2f_ds = downsample_data(Mpr,'time',tsub);
            nnY_ds = quantile(Y_ds(:),0.0005);
            mmY_ds = quantile(Y_ds(:),0.9995);
            nnYf_ds = quantile(Yf_ds(:),0.0005);
            mmYf_ds = quantile(Yf_ds(:),0.99995);
            %%  

            make_avi = true; % save a movie
            if make_avi
                outP = [folder(1:end-10) 'Stablized_AVI/Normcorred_' num2str(clipNum(i)) '_' num2str(iter) '.avi'];
                checkP(outP)
                vidObj = VideoWriter(outP);
                set(vidObj,'FrameRate',30);
                open(vidObj);
            end
            fig = figure;
                screensize = get(0,'Screensize' );
                fac = min(min((screensize(3:4)-100)./[3*d2,d1]),10);
                set(gcf, 'PaperUnits', 'points', 'Units', 'points');
                set(gcf, 'Position', round([100 100 fac*2*d2 fac*d1]));

            for t = 1:1:size(Y_ds,3)
    %             if (0)
                    % plot filtered data
    %                 subplot(131);imagesc(Y_ds(:,:,t),[nnY_ds,mmY_ds]); xlabel('Raw data (downsampled)','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    %                 colormap('bone');
    %                 set(gca,'XTick',[],'YTick',[]);
    %                 subplot(132);imagesc(M1_ds(:,:,t),[nnY_ds,mmY_ds]); xlabel('rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    %                 title(sprintf('Frame %i out of %i',t,size(Y_ds,3)),'fontweight','bold','fontsize',14); 
    %                 colormap('bone')
    %                 set(gca,'XTick',[],'YTick',[]);
% %                     subplot(133);imagesc(M1_ds(:,:,t),[nnY_ds,mmY_ds]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
% %                     colormap('bone')
% %                     set(gca,'XTick',[],'YTick',[]);
    %             else
    %                 % plot full data
                    subplot(121);imagesc(Yf_ds(:,:,t),[nnYf_ds,mmYf_ds]); xlabel('Raw data (downsampled)','fontsize',14,'fontweight','bold'); axis equal; axis tight;
                    colormap('bone');
                    set(gca,'XTick',[],'YTick',[]);
                    subplot(122);imagesc(M1f_ds(:,:,t),[nnYf_ds,mmYf_ds]); xlabel('rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
                    title(sprintf('Frame %i out of %i',t,size(Y_ds,3)),'fontweight','bold','fontsize',14); 
                    colormap('bone')
                    set(gca,'XTick',[],'YTick',[]);
        %             subplot(133);imagesc(M2f_ds(:,:,t),[nnYf_ds,mmYf_ds]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
        %             colormap('bone')
        %             set(gca,'XTick',[],'YTick',[]);
    %             end
                drawnow;
                if make_avi  
                    currFrame = getframe(fig);
                    writeVideo(vidObj,currFrame);    
                end
            end
            if make_avi
                close(vidObj);
                close all
            end
            
            Mr = double(Mr);
            %             if iter==1
            %                 Mr = uint8(Mr);
            outP = [folder(1:end-10) 'Stablized_MAT/Normcorred_' num2str(clipNum(i)) '.mat'];
            checkP(outP)            
            save(outP,'Mr','-v7.3')
%             end
        end
    end
%     MatCat(folder(1:end-10) '\Stablized_MAT');
end
