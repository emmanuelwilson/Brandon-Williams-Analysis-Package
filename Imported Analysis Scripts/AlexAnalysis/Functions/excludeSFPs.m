function excludeSFPs(paths,varargin)
    clc
    close all
    drawnow
    
    warning off all
    if isempty(gcp)
        parpool('local',7);
    end
    pctRunOnAll warning off all
    
    
    fprintf(['Exclude cells based on Spatial Footprints:\n']);
    tmp = [repmat({'\n\t'},[1 length(varargin)]); varargin];
    fprintf(cat(2,tmp{:},'\n\n'))
    
    for p = paths'
        s = load(p{1});
        if isfield(s.processed,'exclude')
%             continue
        end
        
        fprintf(['\t' num2str(p{1}) '\n'])    
        
        SFPs = s.calcium.SFPs;
        SFPs = SFPs./repmat(nanmax(nanmax(SFPs,[],1),[],2),size(SFPs(:,:,1)));
        figure(1)
        imagesc(nanmax(SFPs.^2,[],3))
        BW = roipoly;
        
        
        SFPs = SFPs./repmat(nansum(nansum(SFPs,1),2),[size(SFPs(:,:,1))]);
        
        [x y] = meshgrid(1:size(BW,2),1:size(BW,1));
        
        cent = [permute(nansum(nansum(bsxfun(@times,SFPs,y),1),2),[3 1 2]) ...
            permute(nansum(nansum(bsxfun(@times,SFPs,x),1),2),[3 1 2])];
        
        inds = sub2ind(size(BW),floor(cent(:,1)), ...
            floor(cent(:,2)));
        s.processed.exclude.SFPs = BW(inds);
        fclose all
        save(p{1},'-struct','s','-v7.3');
    end
end
