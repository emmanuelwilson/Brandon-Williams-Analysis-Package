function [sfit sintercept lfit] = fitSigLin(data)
    data(1,:) = fliplr(data(1,:));
    
    sfit = nan(length(data(:,1)),5);
    lfit = nan(length(data(:,1)),1);
    
    x = 1:length(data(1,:));
    for i = 1:length(data(:,1))
        d = cellfun(@nanmedian,data(i,:));
        
        % Sigmoidal fit and plot
        
        [params fval] = fmincon(@(p)sig_error(d,p),[(x(end)./2)+0.5 1 0 2] ,[],[],[],[],...
            [-inf 1e-10 -1 1e-10],[inf 10 1 inf],[],optimoptions('fmincon','Display','none'));
        sfit(i,:) = [fval params];
%         sfit(2,i) = params(1);
        
        hold on
        if i == 1
            plot(fliplr(x),doSig(x,params),'color','r','linestyle',':')
        else
            plot(x,doSig(x,params),'color','b','linestyle',':')
        end
        
%         % Linear fit and plot
%         params = polyfit(x,d,1);
%         pred = polyval(params,x);
% %         if i == 1
% %             plot(fliplr(x),pred,'color','r','linestyle','--')
% %         else
% %             plot(x,pred,'color','b','linestyle','--')
% %         end
%         lfit(i,1) = nanmean((d-pred).^2);
    end
    tmpA = doSig(fliplr(1:0.001:length(data(1,:))),sfit(1,2:5));
    tmpB = doSig(1:0.001:length(data(1,:)),sfit(2,2:5));
    
    tmpFun = @(in,sfit)(doSig(7-in,sfit(1,2:5))-doSig(in,sfit(2,2:5))).^2;
    sintercept = fmincon(@(val)tmpFun(val,sfit),3.5,[],[],[],[],[],[],[],optimoptions('fmincon','Display','none'));
end

function error = sig_error(data,p)
    pred = doSig(1:length(data(1,:)),p);
    error = nanmean((data-pred).^2);
end

function pred = doSig(x,p)
    pred = ([1./(exp((-x+p(1)).*p(2))+1)]./p(4))+p(3);
end