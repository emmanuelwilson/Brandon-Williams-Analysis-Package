function [msout,frameMapout,HeadDirout, HDdegout] = msMerger()
%% Merge ms files with persistent cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function will combine ms files and concactenate them together. The  %
%goal is to stitch together same cells for a prolonged analysis.          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author Emmanuel Wilson

folder = dir('*.mat');     %Looks at all ".mat" files in current directory
frameMapFinal = [];        %Final frame Map variable storage
headTrackFinal = [];       %Final Head coordinate varialbe storage
frameMapbig = [];
HDFinal = [];              %Final Head direction orientation variable strorage
indf = 0;
%Merging all currently useful parts of ms.mat for furthur use
finalStruct = struct();                     %ms.mat final merged version
finalStruct.numFrames = 0;                  %Total number of frames across all videos
finalStruct.mask = [];                      %Mask, the total analysis view within the field of view(FOV) defined by the user in msRun
finalStruct.alignedHeight = 0;              %Aligned height of videos
finalStruct.alignedWidth = 0;               %Aligned width of videos
finalStruct.selectAlignment = 1;            %Alignment selection, user defined in msRun default 1
finalStruct.trace = [];                     %Calcium trace 
finalStruct.firing = [];                    %Firing instances
finalStruct.frameMax = [];                  %Max fluorescence across sessions saved in a single frame
finalStruct.segments = [];                  %cell mask/location in the FOV
finalStruct.frameNum = [];                  %frame number according to video frame was taken from
finalStruct.vidNum = [];                    %Video number/index
finalStruct.vidObj = [];                    %Video itself
finalStruct.hShift = [];                    %height shift/movement correction 
finalStruct.wShift = [];                    %width shift/movement correction
finalStruct.minFrame = {};                  %background frame

% totalFrames = 0;                        %Total number of frames observed (including previous video itterations)
% TimestampMin = 10000000;                %Arbitrary VERY large number for min value storage, must be larger than datenum stamp in folder
% Timestamp = 0;                          %TimeStamp record for indexing
% ind = zeros(length(folder));            %array of indices in chronoligical order
cellIndex = [];                         %Matrix of cell index location across files
count = 1;                              %ms.mat file count
 
for i = 1: length(folder)               %Going through every object in the folder to extract index maping 
    if(contains(folder(i).name,'cellRegistered'))
        cellIndex = load(folder(i).name,'cell_registered_struct');          %extracting the cell index matrix across sessions
        cellIndex = cellIndex.cell_registered_struct.cell_to_index_map;     %Saving the matrix within the function
    end
end

%Merge ms.mat files into one master ms.mat called msfinal.mat
for i = 1 : length(folder)
% %   load(folder(ind(i)).name, 'ms','frameMap','position_track','HD_180_deg_total','SINKdata','HDdeg')        
    load(folder(i).name, 'ms','frameMap','position_track','HD_180_deg_total','SINKdata','HDdeg')
    if (exist('ms','var'))
        finalStruct.numFrames = finalStruct.numFrames + ms.numFrames;
        finalStruct.alignedHeight = ms.alignedHeight;
        finalStruct.alignedWidth = ms.alignedWidth;
        if(count == 1) %initial itteration
            finalStruct.trace = ms.trace;
            finalStruct.firing = ms.firing;
            cellnumber = length(cellIndex(:,1))-length(finalStruct.trace(1,:));
            z = zeros(length(finalStruct.trace(:,1)),cellnumber);
            z3 = zeros(length(ms.segments(:,1,1)),length(ms.segments(1,:,1)),cellnumber);
            finalStruct.trace = [finalStruct.trace z];
            finalStruct.firing = [finalStruct.firing z];
            vidNum = ms.vidNum;
            finalStruct.segments = cat(3, ms.segments,z3) ;
            finalStruct.mask = ms.mask;
            finalStruct.minFrame{1} = ms.minFrame{1};
        else
            vidNum = length(finalStruct.vidNum)+ms.vidNum;
            z = zeros(length(ms.trace(:,1)),length(finalStruct.trace(1,:)));
            finalStruct.trace = [finalStruct.trace; z];
            finalStruct.firing = [finalStruct.firing; z];
            
            if(size(find(finalStruct.mask))< size(find(ms.mask)))
                finalStruct.mask = ms.mask;
            end
            if(sum(sum(finalStruct.minFrame{1})) > sum(sum(ms.minFrame{1})))
                finalStruct.minFrame = ms.minFrame;
            end

        end        
        %cell for cell concactenation
        for j = 1: length(cellIndex(:,1))
            if(cellIndex(j,count) == 0)
                %                     msfinal.trace(length(msfinal.trace(:,j)):(length(msfinal.trace(:,j))+ length(ms.trace(:,j))),j) = 0;
                %                     msfinal.firing(length(msfinal.firing(:,j)):(length(msfinal.firing(:,j))+length(ms.firing(:,j))),j) = 0;
            else
                finalStruct.trace(length(finalStruct.trace(:,1))+1-length(ms.trace(:,1)):length(finalStruct.trace(:,1)),j) = ms.trace(:,cellIndex(j,count));
                finalStruct.firing(length(finalStruct.firing(:,1))+1-length(ms.firing(:,1)):length(finalStruct.firing(:,1)),j) = ms.firing(:,cellIndex(j,count));
                %                msfinal.segments
            end
        end        
        %matrix concactenation
        finalStruct.frameNum = [finalStruct.frameNum ms.frameNum];
        finalStruct.vidNum = [finalStruct.vidNum vidNum];
        finalStruct.vidObj = [finalStruct.vidObj ms.vidObj];
        finalStruct.hShift = [finalStruct.hShift; ms.hShift];
        finalStruct.wShift = [finalStruct.wShift; ms.wShift];
        %         if(isempty(msfinal.minFrame))
        %             msfinal.minFrame = ms.minFrame;
        %         elseif(sum(msfinal.minFrame{1}) > sum(ms.minFrame{1}))
        %             msfinal.minFrame = ms.minFrame;
        %         end
        clear ms
        count = count + 1;
    elseif(exist('frameMap','var'))        
        if isempty(frameMapFinal)
            frameMapFinal = frameMap;
            indf = 1;
        elseif (str2num(folder(i).name(9:end-4)) - indf)>1
            if isempty(frameMapbig)
                frameMapbig = frameMap;
            else
                frameMapbig = [frameMapbig ; max(frameMapbig(:,1))+frameMap];
            end
        else
            frameMapFinal = [frameMapFinal ; max(frameMapFinal(:,1))+frameMap];
            indf = str2num(folder(i).name(9:end-4));
        end
        clear frameMap
    elseif(exist('position_track','var'))
        headTrackFinal = [headTrackFinal; position_track];
        clear position_track
    elseif(exist('SINKdata','var'))
        headTrackFinal = [headTrackFinal; SINKdata];
        clear SINKdata
    end
    if(exist('HDdeg','var'))
        if(isempty(HDFinal))
            HDFinal = HDdeg;
        else        
            HDFinal = [HDFinal; HDdeg];
        end
        clear HDdeg
    elseif(exist('HD_180_deg_total','var'))
        HD_180_deg_total = HD_180_deg_total+180;
        HDFinal = [HDFinal HD_180_deg_total];
        clear HD_180_deg_total
    end
    
end

frameMapFinal = [frameMapFinal ; max(frameMapFinal(:,1))+frameMapbig];
msout = finalStruct;
frameMapout = frameMapFinal;
HeadDirout = headTrackFinal;
HDdegout = HDFinal;

