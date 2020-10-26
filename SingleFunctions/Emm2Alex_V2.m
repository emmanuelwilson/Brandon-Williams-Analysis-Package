%% Emm2Alex will convert all calcium pipeline "ms" files into Alex's "calcium" format with post-processing and SPFs exclusion
%INPUT: -endPath: Will save results in designated folder with the original
%           file folder location as its new file name.
%       -startPath: location of folder with subfolders containing the
%           wanted files to convert


function [] = Emm2Alex_V2(endPath,startPath)

paths = genpath(startPath);
folders = strsplit(paths,';')';
endvars = dir([endPath '/*.mat']);

for i = 20 : length(folders)
    if ~isempty(folders{i})        
        fdate = strsplit(folders{i},'\')';
        fdate = fdate{end-2};
        notprocessed = 0;
        for j = 1 : length(endvars)
            if isempty(find(strcmp(fdate,endvars(j).name),1))
                notprocessed = 1;
                break
            end
        end
        d = dir([folders{i} '/*.mat']);      
        fnames = {d.name};
        if ~isempty(find(strcmp(fnames,'ms.mat'),1)) && notprocessed %&& ~isempty(find(strncmp(fnames,'timestamp',1),1))
            load([folders{i} '/ms.mat']);
            if isfield(ms, 'FiltTraces')
                try
                    if ~isfield(ms , 'deconvolvedSig')
                        display('Deconvolution')
                        if mode(diff(ms.FiltTraces(:,1))) == 0
                            x = [1:length(ms.FiltTraces(:,1))];
                            xq = [1:5:length(ms.FiltTraces(:,1))];
                            temp1 = interp1(x,ms.FiltTraces,xq);
                            ms.FiltTraces = interp1(xq,temp1,x);
                            ms.FiltTraces(find(isnan(ms.FiltTraces))) = 0;
                            save([folders{i} '/ms.mat'],'ms','-v7.3');
                        end
                        ms = msdeconvolve(ms);
                        save([folders{i},'/ms.mat'],'ms','-v7.3');
                    end
                    s.calcium = ms;
                    s.calcium.trace = ms.FiltTraces';
                    s = alignTraceData_noTracking_V2(folders{i},s);
                    s = Emm_converter_excludeSFPs(folders(i),s);
                    s.processed.trace = ms.deconvolvedSig';
                    if ispc
                        varnames = strsplit(folders{i},'\');
                    else
                        varnames = strsplit(folders{i},'/');
                    end
                    if contains(varnames{end},'Miniscope_2')
                        save([endPath,'/',varnames{end-4},'_',varnames{end-3},'_',varnames{end-2},'.mat'],'-struct','s','-v7.3');
                    else
                        save([endPath,'/',varnames{end},'.mat'],'-struct','s','-v7.3');
                    end
                catch
                    fprintf([folders{i}, ' has error process has been terminated and skipping to next session'] )
                end
            end
        end
    end
end