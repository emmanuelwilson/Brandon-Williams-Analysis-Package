%%This Script will convert all of your "msCam.AVI" videos in your directory
%%into one large tiff file.

numFiles = 0;
o = [];
filePrefix = 'behavCam';
% find avi files
aviFiles = dir([pwd '\*.avi']);

%find the total number of relevant video files
for i=1:length(aviFiles)
    endIndex = strfind(aviFiles(i).name,'.avi');        %find the name of current .avi file
    if (~isempty(strfind(aviFiles(i).name,filePrefix)))
        numFiles = max([numFiles str2double(aviFiles(i).name((length(filePrefix)+1):endIndex))]);     % +1 count for relevant .avi files
        %         o(length(o)+1) = i;
    end
end
Timestampmin = 1000000;     %Arbitrary VERY large number for min value storage, must be larger than datenum stamp in folder
Timestamp = 0;              %Total number of frames observed (including previous video itterations)
o = NaN(1,numFiles);     %index movie file order
anomilynames = NaN(1,4);
count = 0;

%Sorting through file timestamps and marking the index's for future reference
for i = 1 : numFiles
    for j = 1 : length(aviFiles)
        if(aviFiles(j).datenum == Timestamp)
            count = count +1;
            anomilynames(1,count) = j;
        end
        if (~isempty(strfind(aviFiles(j).name,filePrefix)) && aviFiles(j).datenum < Timestampmin && aviFiles(j).datenum > Timestamp)
            o(1,i) = j;                         %Store file location within folder
            Timestampmin = aviFiles(j).datenum;   %Set minimum time to lowest value found thought the loop
            %                 prevname = folder(j).name;
        end
        if count > 2
            if aviFiles(anomilynames(1,2)).bytes < aviFiles(anomilynames(1,3)).bytes
                o(1,i) = anomilynames(1,3);                         %Store file location within folder
                Timestampmin = aviFiles(j).datenum;   %Set minimum time to lowest value found thought the loop
            end
        end
    end
    Timestamp = Timestampmin;                   %Reset loop boundaries
    Timestampmin = Timestampmin*100;             %Reset loop boundaries
    count = 0;
    
end
if length(o)>1 && aviFiles(o(1,length(o(1,:)))).bytes > aviFiles(o(1,length(o(1,:))-1)).bytes
    temp = o(1,length(o(1,:)));
    o(1,length(o(1,:)))= o(1,length(o(1,:))-1);
    o(1,length(o(1,:))-1)= temp;
end
for i =1 : numFiles
    obj = VideoReader(aviFiles(o(i)).name);
    vid = read(obj);
    frames = obj.NumberOfFrames;    
    for x = 1 : frames
        if i == 1 && x ==1
            imwrite(vid(:,:,:,x),'msCam.tif');
        else
            imwrite(vid(:,:,:,x),'msCam.tif','WriteMode','append');
        end
    end
end
