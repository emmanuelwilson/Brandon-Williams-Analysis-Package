function [vecs pvals] = getBurstVectors(p,T,pause_thresh,cell_thresh,isPC)
    vecs = [];
    pvals = [];

    minShift = 900;
    sampling = 30;

    v = [0 sqrt(nansum(diff(p,[],2).^2,1))].*sampling;
%     sp = imfilter(p,fspecial('gauss',[1 5.*sampling],sampling.*0.2),'same','replicate');
%     sv = [0 sqrt(nansum(diff(sp,[],2).^2,1))].*sampling;
%     v = sv;
% % %     
% % %     nsim = 100;
% % %     null = nan(length(T(:,1))+1,nsim);
% % %     for si = 1:nsim
% % %         gT = circshift(T,[0 minShift+randi(length(T(1,:))-minShift.*2)]);
% % %         null(:,si) = cumsum(histc(nansum(gT(:,v<pause_thresh)>0,1),[0:1:length(T(:,1))]));
% % % 
% % % %         cumHist([{nansum(T(:,v<pause_thresh)>0,1)}; {nansum(gT(:,v<pause_thresh)>0,1)}],[0:1:60])
% % % %         drawnow
% % %     end
% % %     null = null./null(end,end);
% % %     actual = cumsum(histc(nansum(T(:,v<pause_thresh)>0,1),[0:1:length(T(:,1))])');
% % %     actual = actual./actual(end,end);
% % %     
% % %     figure
% % %     set(gcf,'position',[50 50 1200 800])
% % %     subplot(1,2,1)
% % %     plot(actual)
% % %     hold on
% % %     plot(nanmedian(null,2))
% % %     subplot(1,2,2)
% % %     imagesc(actual>=nanmedian(null,2))
% % %     drawnow

%     tmp = nanmean(T(isPC,:)>0,1);
%     figure(1)
%     mkGraph([{tmp(v<3)'} {tmp(v>=3)'}])
%     drawnow
    
    cT = imfilter(T,fspecial('gauss',[1 300],30),'same','replicate');
    vecs = cT(:,nanmean(cT>0.001,1)>=cell_thresh & v<pause_thresh);
    
end