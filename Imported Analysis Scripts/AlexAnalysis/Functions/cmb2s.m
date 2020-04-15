function cmb2s(inputFolder,outputFolder,overwrite)

    for i = 1:length(objectDirs)
        isAnimal = ismember(cat(1,objectDirs{:,1}),objectDirs{i,1});
        isAlreadyLoaded = [outputFolder '/' num2str(objectDirs{i,1}) '/' ...
                objectDirs{i,5}(find(ismember(objectDirs{i,5},'\'),1,'last')+1:end) '/data.mat'];
        if exist(isAlreadyLoaded,'file') == 2 && ~overwrite% If struct already exists, skip it!
            continue
        end

        s = struct;

        objString = [strcat(objectDirs{i,5},'\batch_object.mat')];
        tObjString = objString;
        tObjString(ismember(objString,'\')) = '/';
        fprintf(['\n\t\tConverting Objects to Structs:  ' tObjString ' (' num2str(i) ' of ' num2str(length(objectDirs))  ')'])

        if ismember(i,[918 1060 1129]) %%% Catch for random errors. 366 Timestamps not monotonically increasing
            continue
        end

        objString = [strcat(objectDirs{i,5},'/batch_object.mat')];
        if exist(objString,'file') ~= 2
            if exist([objectDirs{i,5} '\ClusterCut'],'dir') == 7
                files = dir([objectDirs{i,5} '\ClusterCut']);
                if length(files) > 2

                    curDir = pwd;
                    MakeObject_batch(objectDirs{i,5},false)
                    cd(curDir);
                else
                    current = pwd;
                    cd(objectDirs{i, 5})

                    if exist('TT1.mat') == 2, TT1 = load('TT1.mat'); Sc1 = TT1(:,1:2); cd(strcat(objectDirs{i, 5},'\ClusterCut')),save Sc1 Sc1,cd(objectDirs{i, 5});end
                    if exist('TT2.mat') == 2, TT2 = load('TT2.mat'); Sc2 = TT2(:,1:2); cd(strcat(objectDirs{i, 5},'\ClusterCut')),save Sc2 Sc2,cd(objectDirs{i, 5});end 
                    if exist('TT3.mat') == 2, TT3 = load('TT3.mat'); Sc3 = TT3(:,1:2); cd(strcat(objectDirs{i, 5},'\ClusterCut')),save Sc3 Sc3,cd(objectDirs{i, 5});end
                    if exist('TT4.mat') == 2, TT4 = load('TT4.mat'); Sc4 = TT4(:,1:2); cd(strcat(objectDirs{i, 5},'\ClusterCut')),save Sc4 Sc4,cd(objectDirs{i, 5});end

                    if exist('TT1.txt') == 2, TT1 = load('TT1.txt'); Sc1 = TT1(:,1:2); cd(strcat(objectDirs{i, 5},'\ClusterCut')),save Sc1 Sc1,cd(objectDirs{i, 5});end
                    if exist('TT2.txt') == 2, TT2 = load('TT2.txt'); Sc2 = TT2(:,1:2); cd(strcat(objectDirs{i, 5},'\ClusterCut')),save Sc2 Sc2,cd(objectDirs{i, 5});end 
                    if exist('TT3.txt') == 2, TT3 = load('TT3.txt'); Sc3 = TT3(:,1:2); cd(strcat(objectDirs{i, 5},'\ClusterCut')),save Sc3 Sc3,cd(objectDirs{i, 5});end
                    if exist('TT4.txt') == 2, TT4 = load('TT4.txt'); Sc4 = TT4(:,1:2); cd(strcat(objectDirs{i, 5},'\ClusterCut')),save Sc4 Sc4,cd(objectDirs{i, 5});end

                    files = dir([objectDirs{i,5} '\ClusterCut']);
                    if length(files) <= 2
                        missingCuts = [missingCuts; objectDirs(i,5)];
                        continue
                    end
                    MakeObject_batch(objectDirs{i,5},false)
                    cd(current)

                    curDir = pwd;
                    MakeObject_batch(objectDirs{i,5},false)
                    cd(curDir);
                end
            else
                mkdir(strcat(objectDirs{i,5},'\ClusterCut'))
                current = pwd;
                cd(objectDirs{i, 5})

                if exist('TT1.mat') == 2, TT1 = load('TT1.mat'); Sc1 = TT1(:,1:2); cd(strcat(objectDirs{i, 5},'\ClusterCut')),save Sc1 Sc1,cd(objectDirs{i, 5});end
                if exist('TT2.mat') == 2, TT2 = load('TT2.mat'); Sc2 = TT2(:,1:2); cd(strcat(objectDirs{i, 5},'\ClusterCut')),save Sc2 Sc2,cd(objectDirs{i, 5});end 
                if exist('TT3.mat') == 2, TT3 = load('TT3.mat'); Sc3 = TT3(:,1:2); cd(strcat(objectDirs{i, 5},'\ClusterCut')),save Sc3 Sc3,cd(objectDirs{i, 5});end
                if exist('TT4.mat') == 2, TT4 = load('TT4.mat'); Sc4 = TT4(:,1:2); cd(strcat(objectDirs{i, 5},'\ClusterCut')),save Sc4 Sc4,cd(objectDirs{i, 5});end

                if exist('TT1.txt') == 2, TT1 = load('TT1.txt'); Sc1 = TT1(:,1:2); cd(strcat(objectDirs{i, 5},'\ClusterCut')),save Sc1 Sc1,cd(objectDirs{i, 5});end
                if exist('TT2.txt') == 2, TT2 = load('TT2.txt'); Sc2 = TT2(:,1:2); cd(strcat(objectDirs{i, 5},'\ClusterCut')),save Sc2 Sc2,cd(objectDirs{i, 5});end 
                if exist('TT3.txt') == 2, TT3 = load('TT3.txt'); Sc3 = TT3(:,1:2); cd(strcat(objectDirs{i, 5},'\ClusterCut')),save Sc3 Sc3,cd(objectDirs{i, 5});end
                if exist('TT4.txt') == 2, TT4 = load('TT4.txt'); Sc4 = TT4(:,1:2); cd(strcat(objectDirs{i, 5},'\ClusterCut')),save Sc4 Sc4,cd(objectDirs{i, 5});end

                files = dir([objectDirs{i,5} '\ClusterCut']);
                if length(files) <= 2
                    missingCuts = [missingCuts; objectDirs(i,5)];
                    continue
                end

                cd(current)

                curDir = pwd;
                MakeObject_batch(objectDirs{i,5},false)
                cd(curDir);
            end
        end
        if exist(objString,'file') == 2
            load(objString);
            isAnimal = ismember(cat(1,objectDirs{:,1}),objectDirs{i,1});

            s.properties.rat = num2str(objectDirs{i,1});
            s.properties.type = objectDirs{i,2};
            s.properties.age = objectDirs{i,3};
            s.properties.environment = objectDirs{i,4};
            s.properties.date = objectDirs{i,5}(find(ismember(objectDirs{i,5},'\'),1,'last')+1:end);
            s.properties.sessionNum = objectDirs{i,6}; %find(ismember(objectDirs(isAnimal,5),objectDirs(i,5)));


            s.pos.p = [root.x'; root.y'] .* 0.21;
            s.pos.ts = root.ts';
            s.pos.hd = deg2rad(root.headdir');
            s.pos.vel = root.vel';

            s.p.cm2pix = 2.5;
            s.p.kern_sd = 4;

            s.unit = repmat({[]},size(root.spike));
            for t = 1:length(s.unit(:,1))
                for c = 1:length(s.unit(1,:))
                    if ~isempty(root.spike(t,c).ts)
                        s.unit{t,c} = root.spike(t,c).ts; 
                    end
                end
            end

%             for channel = 1:length(root.b_lfp)
%                 s.lfp(channel).signal = root.b_lfp(channel).signal;
%                 s.lfp(channel).theta = root.b_lfp(channel).theta;
%                 s.lfp(channel).theta_phase = root.b_lfp(channel).theta_phase;
%                 s.lfp(channel).theta_amplitude = root.b_lfp(channel).theta_amplitude;
%             end

            outP = [outputFolder '/' s.properties.rat '/' objectDirs{i,5}(find(ismember(objectDirs{i,5},'\'),1,'last')+1:end) '/data.mat'];
            checkP(outP);
            save(outP,'-struct','s','-v7.3')
        else
            
        end
    end

end