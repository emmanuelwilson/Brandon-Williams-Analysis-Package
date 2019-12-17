%% Crawl script for deconvolved ms structure creation

clear all

folder = dir(pwd);
oldCD = pwd;
for i = 3 : length(folder)
    if folder(i).isdir == 1
        subdir = [pwd,'\',folder(i).name];
        subfolder = dir(subdir);
        fnames = {subfolder.name};
        if ~isempty(find(strcmp(fnames,'ms.mat'),1)) && ~isempty(find(strcmp(fnames,'SFP.mat'),1)) && ~isempty(find(strncmp(fnames,'timestamp.dat',1),1)) && isempty(find(strcmp(fnames,'msDeconvolved.mat'),1))
            cd([pwd,'/',folder(i).name]);
            load('ms.mat')
            x = [1:length(ms.FiltTraces(:,1))];
            xq = [1:5:length(ms.FiltTraces(:,1))];
            temp1 = interp1(x,ms.FiltTraces,xq);
            ms.FiltTraces = interp1(xq,temp1,x);
            ms.FiltTraces(find(isnan(ms.FiltTraces))) = 0;
            ms = msdeconvolve(ms);
            save('msDeconvolved','ms')
            cd(oldCD);
        end
    end
end