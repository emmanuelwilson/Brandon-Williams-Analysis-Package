function normTrace(paths)
    clc
    fprintf(['\t\nNorming traces with 2nd-order autoregressive model.\n\n'])
    warning off all
    for p = paths'
        s = load(p{1});
%         s.processed.norming.trace = s.processed.trace;
%         s.processed.norming.traceModel = s.processed.traceModel;
        fprintf(['\t\t' num2str(p{1}) '\n\t\t\tModeling spike trains:      '])
        
        t = s.calcium.FiltTraces';
        isGood = find(diff(t(1,:),[],2)~=0);
        for k = 1:length(t(:,1))
            t(k,isGood(1):isGood(end)) = linterp(isGood,t(k,isGood),isGood(1):isGood(end));
        end
        
%         zt = bsxfun(@rdivide,bsxfun(@minus,t,nanmean(t,2)),nanstd(t,[],2));
%         dzt = [nan(length(zt(:,1)),1) diff(zt,[],2)];
%         isRise = dzt > 0 & zt > 1;
% %         isRise = [zeros(length(isRise(:,1)),1) diff(isRise,[],2)==1];
%         s.processed.trace = isRise;
        
         t = s.calcium.FiltTraces';
%         t = s.calcium.FiltTraces(s.processed.validTraceFrames(:,1),:)';
        isGood = find(diff(t(1,:),[],2)~=0);
        for k = 1:length(t(:,1))
            t(k,isGood(1):isGood(end)) = linterp(isGood,t(k,isGood),isGood(1):isGood(end));
        end
        
        isRise = nan(size(t));
        t = detrend(t')'; % norm to the median
%         parfor k = 1:length(t(:,1))
%             bsd = nanstd(t(k,t(k,:)<0))./sqrt(1-(2./pi));
%             zt = t(k,:)./bsd;
%             dzt = [nan diff(zt,[],2)];
%             isRise(k,:) = dzt > 0 & zt > 3;
%         end
%         s.processed.trace = isRise;
%         s.processed.traceModel = {'Binarization'};
        
        
%         s.processed.norming(2).trace = s.processed.trace;
%         s.processed.norming(2).traceModel = s.processed.traceModel;
        
        t = t';
        spikes = nan(size(s.calcium.FiltTraces));
        for i = 1:length(s.calcium.FiltTraces(1,:))
            str_f = sprintf('%6.1f',100.*i/length(s.calcium.FiltTraces(1,:)));
            fprintf([repmat('\b',[1 6]) str_f])
            [blah spikes(:,i)] = deconvolveCa(t(:,i),'type','ar2','b','optimize_b');
        end
        fprintf('\n')
        s.processed.trace = spikes';
        s.processed.traceModel = {'AutoReg2'};
%         
        save(p{1},'-struct','s','-v7.3');
    end
end