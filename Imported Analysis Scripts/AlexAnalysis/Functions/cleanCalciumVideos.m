function cleanCalciumVideos(p)
    clc
    
    fprintf('Checking and correcting blank frames... \n')
    
    % make sure there aren't any pure-black frames, which would crash msRun
    files = dir([p '/msCam*.avi']);
    names = {files(:).name};
    for name = names
        obj = VideoReader([p '/' name{1}]);
        frames = [];
        
        tmp = read(obj);
        tmp2 = reshape(tmp,[obj.Width.*obj.Height length(tmp(1,1,1,:))]);
        
        interpFrames = find(nanmean(tmp2)<=5);
        
        if ~isempty(interpFrames)
            fprintf(['\tFrames interpolated for:  ' num2str(length(interpFrames)) '\n' p '/' name{1} '\n'])
            for i = interpFrames
                if i==1
                    tmp(:,:,:,i) = tmp(:,:,:,i+1);
                elseif i == length(tmp(1,1,1,:))
                    tmp(:,:,:,i) = tmp(:,:,:,i-1);
                else
                    tmp(:,:,:,i) = tmp(:,:,:,i-1);
                end
            end


            outObj = VideoWriter([p '/' name{1}],'Grayscale AVI');
            open(outObj)
            writeVideo(outObj,tmp);
            close(outObj)
        end
    end
end