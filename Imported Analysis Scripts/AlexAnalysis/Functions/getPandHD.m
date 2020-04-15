function getPandHD(paths)
    clc
    close all
    drawnow
    
    params = ebcMapParams();
       
    warning off all
    if isempty(gcp)
        parpool('local',7);
    end
    pctRunOnAll warning off all
    
    
    envSizes = [60 60; 75 75];
    fprintf(['Computing P and HD for 2 LEDs:\n']);
    for p = paths'
        s = load(p{1});
        s.processed.envSize = envSizes((s.properties.session(1)=='B')+1,:);
        
        fprintf(['\n\t' num2str(p{1})]) 
        if isfield(s.pos,'mp')
            mp = s.pos.mp;
            tp = mp(:,:,2)-mp(:,:,1);
            [t d] = cart2pol(tp(1,:),tp(2,:));
            t(d<10 | d>25) = nan;
            s.processed.hd = rad2deg(mod(interpCircNaNs(t')'-(pi./2),2.*pi));
            
            tp = nanmean(mp,3);
            tp(:,d<10 | d>25) = nan; %%%%% CHANGES WITH VIDEO SIZE; COULD NORM TO ENV SIZE
            tp = interpNaNs(tp')';
            
            mm = nanmin(nanmin(mp,[],2),[],3);
            mx = nanmax(nanmax(bsxfun(@minus,mp,mm),[],2),[],3);
            
            s.processed.p = bsxfun(@times,bsxfun(@rdivide,bsxfun(@minus,tp,mm),mx),s.processed.envSize');
            
%             md = cart2pol(diff(s.processed.p(1,:)),diff(s.processed.p(2,:)));
%             scatter(md,s.processed.hd(2:end))
%             
%             amt = length(s.processed.hd);
%             figure
%             plot(s.processed.p(1,1:amt),s.processed.p(2,1:amt),'color',[0.5 0.5 0.5])
%             hold on
%             plot([s.processed.p(1,1:amt); s.processed.p(1,1:amt)+sind(s.processed.hd(1,1:amt)).*0.5],...
%                 [s.processed.p(2,1:amt); s.processed.p(2,1:amt)+cosd(s.processed.hd(1,1:amt)).*0.5],...
%                 'color','r')
            
            save(p{1},'-struct','s','-v7.3');
        end
    end
end