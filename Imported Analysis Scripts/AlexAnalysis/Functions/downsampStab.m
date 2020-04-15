function downsampStab(folder)
    files = getFilePaths('StablizedVideos','.avi');
    
    ds = 10;
    for p = files'
        a = read_file(p{1});
        
        outObj = VideoWriter([p{1}(1:end-4) '_ds']);
        open(outObj)
        writeVideo(outObj,permute(a(:,:,1:ds:end),[1 2 4 3]));
        close(outObj)
    end
end