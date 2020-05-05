function mapAndPreprocess(paths)
    clc
    envSizes = [60 60; 75 75];
    for i = 1:length(paths)
        fprintf(['\n\t' paths{i}])
        
        s = load(paths{i});
        
        s.processed.p = s.pos.p(:,s.frameMap);
        s.processed.hd = s.pos.hd(s.frameMap);
        s.processed.envSize = envSizes((s.properties.session(1)=='B')+1,:);
        
        
        s.processed.p = s.processed.envSize(1).*[(s.processed.p - nanmin(s.processed.p(:))) ./ ...
            (nanmax(s.processed.p(:)) -nanmin(s.processed.p(:)))];
        
        s.processed.p = interpNaNs(s.processed.p')';

        MD = cart2pol(diff( s.processed.p(1,:),[],2),diff( s.processed.p(2,:),[],2));
        MD(end+1) = MD(end);

        s.processed.md = rad2deg(MD);
        
%         shift = (s.processed.envSize'./2) - nanmedian(s.processed.p,2);
%         s.processed.p = bsxfun(@minus,s.processed.p,shift); 
        
        s.processed.hd  = interpCircNaNs(s.processed.hd')';
        
        save(paths{i},'-struct','s','-v7.3')
    end
end