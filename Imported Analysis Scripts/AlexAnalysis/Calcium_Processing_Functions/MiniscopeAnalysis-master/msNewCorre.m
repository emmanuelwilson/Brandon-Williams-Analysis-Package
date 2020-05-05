function ms = msNewCorre(ms,isnonrigid,doWindowed,varargin)

warning off all

    separator = '\'; % For pc operating systems

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
            nanmax(ceil(limit(1)),1):nanmin(ceil(limit(1)+limit(3)),length(Yf(:,1,1))),:);
    end
    bound = round(nanmin(size(Yf(:,:,1)))/2);
    Y = imfilter(Yf,psf,'symmetric');
    [d1,d2,T] = size(Y);
    
    if isempty(mask)
        figure(1)
        imagesc(Y(:,:,1));
        colormap parula
        colorbar
        mask = getrect();
    end
    
    fftalign(a,b,false)
       
end

close(writerObj);

end