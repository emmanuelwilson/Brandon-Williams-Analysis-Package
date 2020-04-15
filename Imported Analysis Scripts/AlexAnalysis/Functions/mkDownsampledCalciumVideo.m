function mkDownsampledCalciumVideo(folder)
%     folder = 'Data/OpenFieldData/ANP021117/18-04-26_Baseline/19-May-2018_8-42-15/msvideo.avi';
%     folder = 'X:/Alex/ThreeRooms_TwoDoors_OneExperiment/Data/RawCalcium_1LED/AKCA106/18-05-13_Baseline';
%     
%     files = dir([folder '/msCam*.avi']);
%     
%     clipNum = nan(1,length(files));
%     for i = 1:length(files)
%         clipNum(i) = str2num(files(i).name(6:end-4));
%     end
%     [a b] = sort(clipNum);
%     files = files(b);
%     
%     names = {files(:).name};
%     ds = 4;
%     allF = [];
%     for name = names(3:7)
%         obj = VideoReader([folder '/' name{1}]);
%         tmp = read(obj);
%         new = uint8(zeros([round(size(tmp(:,:,1,1))./2) round(length(tmp(1,1,1,:))./ds)]));
%         for i = 1:ds:length(tmp(1,1,1,:))
%             new(:,:,(i-1)./ds + 1) = imresize(tmp(:,:,1,i),[round(size(tmp(:,:,1,1))./2)]);
%         end    
%         allF = cat(3,allF,new);
%     end
%     
%     outObj = VideoWriter(['AKCA106_ForShow.avi'],'Grayscale AVI');
%     open(outObj)
%     writeVideo(outObj,allF(:,:,1:1000));
%     close(outObj)
%     
    
%     backgroundImage = VideoReader('./Data/RawCalcium_1LED/ANP030803/18-04-28_Baseline/msCam10.avi');
%     backgroundFrame = read(backgroundImage,500);
%     

%%%%%%%%%%%%%%%%%%

    file = 'msvideo.avi';
    obj = VideoReader(file);
    allF = read(obj);
    fclose all;

    
    s = load('./MatlabData_OpenFieldData/OpenFieldData/ANP030803/18-04-28_Baseline');
    tmp = s.calcium.SFPs;
    
    
%     numColors = 32;
%     doColors = hsv(numColors);
%     forPrettyPlot = zeros([size(tmp(:,:,1)) 3]);
%     for i = 1:length(tmp(1,1,:))
%         tmp2 = tmp(:,:,i);
%         tmp2 = tmp2./nanmax(tmp2(:));
% %         tmp2 = tmp2.^2;
%         tmp2(tmp2<0.5) = 0;
%         forPrettyPlot = forPrettyPlot + bsxfun(@times,tmp2,permute(doColors(mod(i,numColors-1)+1,:),[3 1 2]));
%     end
%     figure(1)
%     set(gcf,'position',[50 50 700 300])
%     subplot(1,2,1)
%     file = 'Data/OpenFieldData/ANP030803/18-04-28_Baseline/19-May-2018_10-51-2/msvideo.avi';
%     obj = VideoReader(file);
%     allF = read(obj,100);
%     fclose all;
%     imshow(allF)
%     subplot(1,2,2)
%     imshow(forPrettyPlot)
    
    
%     shift = round([[obj.Height - size(tmp,1)] [obj.Width - size(tmp,2)]]./2);
%     tmp = tmp(shift(1):end-shift(1)+1,shift(2):end-shift(2)+1,:);
    isGood = find(s.processed.splithalf.p <= 0.05);
    choose = isGood([94 23 49 64 37 26 146 1 156 103 193 150 93]); % [14 16 18 19 83]    14 16
    filter = zeros([round(size(allF(:,:,1))) length(choose)]);
    for i = 1:length(choose)
        good = tmp(:,:,choose(i));
        good = good./nanmax(good(:));
%         good(good < 0.5) = 0;
% %         good = imresize(tmp(:,:,choose(i)),size(filter(:,:,1)));
%         good = circshift(good,[2 0]) + circshift(good,[-2 0]) + ...
%             circshift(good,[0 2]) + circshift(good,[0 -2]);
        filter(:,:,i) = good;
    end
    doColors = abs(1-(hsv(length(filter(1,1,:)))));
%     doColors = doColors(randperm(length(doColors)),:);
%     doTrace = s.processed.trace(choose,1:10:end);
    step = 1;
    doTrace = s.calcium.FiltTraces';
%     doTrace = s.calcium.trace';
    doTrace = doTrace(choose,1:step:end);
    doTrace = bsxfun(@minus,doTrace,min(doTrace,[],2));
    doTrace = bsxfun(@rdivide,doTrace,max(doTrace,[],2));
%     doTrace = bsxfun(@minus,doTrace
    
    pretty = zeros([size(filter(:,:,1)) 3]);
    separate = zeros([size(filter(:,:,1)) 3 length(choose)]);
    for i = 1:length(choose)
        pretty = pretty + bsxfun(@times,filter(:,:,i),permute(abs(1-doColors(i,:)),[1 3 2]));
        separate(:,:,:,i) = bsxfun(@times,filter(:,:,i),permute(abs(1-doColors(i,:)),[1 3 2]));
    end
%     pretty = abs(1-pretty);
%     pretty = imfilter(double(pretty),fspecial('gauss',[31 31],2),'same');
    
    
%     new = zeros([size(filter(:,:,1)) 3 length(doTrace(1,:))]);
%     for i = 1:length(doTrace(:,1))
% %         new = new + bsxfun(@times,permute(bsxfun(@times,...
% %             filter(:,:,i),permute(doTrace(i,:),[1 3 2])),[1 2 4 3]),pretty);
%         new = new + bsxfun(@times,permute(bsxfun(@times,...
%             filter(:,:,i),permute(doTrace(i,:),[1 3 2])),[1 2 4 3]),repmat(permute(doColors(i,:),[1 3 2]),[size(filter(:,:,1))]));
%     end
%     
%     i = length(new(1,1,1,:))-1;
    
%     allF = double(allF);
%     allF = bsxfun(@minus,allF,nanmin(allF,[],4));
%     allF = bsxfun(@rdivide,allF,nanmax(allF,[],4));

    movObj = getframe();
    count = 0;
%     minF = double(nanmin(allF,[],4));
%     minF = double(nanmean(allF,4));
%     maxF = double(nanmax(allF,[],4));
%     minF = minF+(maxF-minF).*0.2;
    step = 10;
    for i = [1:step:length(allF)] % [1:6:1200 1206:12:1500 1512:60:2400 2460:240:length(allF(1,1,1,:))-1]
        count = count+1;
        figure(1)
        set(gcf,'position',[50 50 750 500],'color','k')
%         subplot(2,1,1)
%         doF = abs((double(allF(:,:,:,i))-minF)./(maxF-minF)); % ./255
%         forImage = (abs(1-pretty).*(repmat(doF,[1 1 3])));
%         forImage = (abs(pretty).*(repmat(doF,[1 1 3])));
%         h = image(abs(1-forImage));
%         axis off 
%         set(h,'alphadata',doF);
%         set(gca,'xtick',[],'ytick',[],'color','k')
%         axis equal
%         drawnow
        

%         subplot(2,2,1)
        imshow(allF(:,:,1,i))
%         imagesc((double(allF(:,:,1,i+1))-double(allF(:,:,1,i))))
%         colormap gray
%         
%         subplot(2,2,3)
%         doF = zeros([size(filter(:,:,1,1)) 3]);
%         for j = 1:length(choose)
%             doF = doF + separate(:,:,:,j).*doTrace(j,i);
%         end
%         image(doF)
% %         set(gca,'ydir','normal')
%         axis equal
%         axis off
%         drawnow
%         
% %         colormap gray
% %         hold on
% %         h = imagesc(new(:,:,:,i));
% %         set(h,'alphadata',(double(~all(new(:,:,:,i)==0,3))));
% %         axis tight
% %         axis off  
% 
%         subplot(1,2,2)
%         rectangle('position',[0 0 40 40],'facecolor','k');
%         hold on
%         plot(s.processed.p(1,1:i.*step+1),...
%             s.processed.p(2,1:i.*step+1),'color',[0.8 0.8 0.8])
%         for j = 1:length(choose)
%             plot(s.processed.p(1,logical(s.processed.trace(choose(j),1:i.*step+1))),...
%                 s.processed.p(2,logical(s.processed.trace(choose(j),1:i.*step+1))),'color',abs(1-doColors(j,:)),...
%                 'markerfacecolor',abs(1-doColors(j,:)),'linestyle','none','marker','o','markersize',3);
%         end
%         hold off
%         axis equal
%         axis off
%         set(gca,'xlim',[0 40],'ylim',[0 40])
%         drawnow
%         
%         pause(0.1)
%         
        movObj(count) = getframe(gcf);
        close all
    end

    movie2avi(movObj(2:end-1),'AKCA3D05_ForShow','fps',30)
    
    outObj = VideoWriter(['AKCA3D05_ForShow.avi'],'Grayscale AVI');
    open(outObj)
    writeVideo(outObj,movObj);
    close(outObj)

    
    
%     s = load('./MatlabData/RawCalcium_1LED/AKCA106/18-05-13_Baseline');
%     tmp = s.calcium.segments;
%     shift = round([[obj.Height - size(tmp,1)] [obj.Width - size(tmp,2)]]./2);
%     tmp = tmp(shift(1):end-shift(1)+1,shift(2):end-shift(2)+1,:);
%     isGood = find(s.processed.splithalf.p <= 0.05);
%     choose = isGood([1:1:40]);
%     filter = zeros(size(allF(:,:,1)));
%     for i = 1:length(choose)
%         good = imresize(tmp(:,:,choose(i)),size(allF(:,:,1)));
%         good = circshift(good,[2 0]) + circshift(good,[-2 0]) + ...
%             circshift(good,[0 2]) + circshift(good,[0 -2]);
%         filter(good~=0) = i;
%     end
%     doColors = hsv(nanmax(filter(:)));
%     pretty = zeros([size(filter) 3]);
%     for i = 1:length(choose)
%         pretty = pretty + bsxfun(@times,filter==i,permute(doColors(i,:),[1 3 2]));
%     end
%     pretty = imfilter(double(pretty),fspecial('gauss',[31 31],1),'same');
%     s.processed.trace(choose,1:10:end);
    
%     colorizedF = double(allF);
%     colorizedF = imfilter(double(colorizedF),fspecial('gauss',[31 31],3),'same');
% %     colorizedF = permute(imfilter(permute(colorizedF,[1 3 2]),fspecial('gauss',[1 31],5),'same'),[1 3 2]);
%     colorizedF = bsxfun(@minus,colorizedF,nanmin(colorizedF,[],3));
%     colorizedF = bsxfun(@rdivide,colorizedF,nanmax(colorizedF,[],3));
%     colorizedF = permute(colorizedF,[1 2 4 3]);
%     colorizedF = bsxfun(@times,colorizedF,pretty);
%     
%     for i = 1:length(colorizedF(1,1,1,:))
%         imshow(colorizedF(:,:,:,i))
%         axis tight
%         axis off
%         drawnow
%     end
%     
%     outObj = VideoWriter(['AKCA106_ForShow.avi'],'Grayscale AVI');
%     open(outObj)
%     writeVideo(outObj,allF(:,:,1:1000));
%     close(outObj)
    
end