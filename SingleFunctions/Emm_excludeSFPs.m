function [] = Emm_excludeSFPs(p,varargin)
clc
close all
drawnow

paths = genpath(p);
if ispc
    folders = strsplit(paths,';')';
else
    folders = strsplit(paths,':')';
end
for i = 1 : length(folders)
    if ~isempty(folders{i})        
        d = dir(folders{i});
        fnames = {d.name};
        if ~isempty(find(strcmp(fnames,'ms.mat'),1))
            cd(folders{i});
            load('ms.mat')
            if ~isfield(ms,'exclude')                                
%                 if isempty(gcp)
%                     parpool('local',7);
%                 end                         
                
                fprintf(['Exclude cells based on Spatial Footprints:\n']);
                                
                fprintf(['\t' ms.dirName '\n'])
                
                SFPs = ms.SFPs;
                SFPs = SFPs./repmat(nanmax(nanmax(SFPs,[],1),[],2),size(SFPs(:,:,1)));
                figure(1)
                imagesc(nanmax(SFPs.^2,[],3))
                BW = roipoly;
                
                SFPs = SFPs./repmat(nansum(nansum(SFPs,1),2),[size(SFPs(:,:,1))]);
                
                [x y] = meshgrid(1:size(BW,2),1:size(BW,1));
                
                cent = [permute(nansum(nansum(bsxfun(@times,SFPs,y),1),2),[3 1 2]) ...
                    permute(nansum(nansum(bsxfun(@times,SFPs,x),1),2),[3 1 2])];
                
                inds = sub2ind(size(BW),floor(cent(:,1)), floor(cent(:,2)));
                ms.exclude.SFPs = BW(inds);
                fclose all
                save('ms','ms','-v7.3');
            end
        end
    end
end
end
