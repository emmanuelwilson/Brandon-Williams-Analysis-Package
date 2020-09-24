function [] = loopGraphs_v3(p, init)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Takes the initial 3 letters of your animal ID as input to find all 'ms.mat' files, 
%Attempt to create and save .fig and .png files of SFP and FiltTraces of recording. 
%Displays failure message otherwise, for each file. 
%All files saved to same directory as ms.mat file.
%
%Example use: loopGraphs_v2(pwd, 'ZHA')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% if you WANT TO SHOW GRAPHS then set 'show' to 'true' %%%
%%% if you don't want to save figures, set 'save' to 'false' %%%


% how long is your animal ID? For example, 'ZHA001' is 6 characters long %

% What file format do you want to save as? .fig is saved
% automatically change '.svg' if you'd prefer a different image format

%lastly, what range should the colorbar display? set to 0.09 to 0.3 by default

show = true;
save = true;

animalID_length = 9;
fileformat = '.svg';

cmin = 0.09;
cmax = 0.3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','matlab:audiovideo:VideoReader:FileNotFound')
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
        if isempty(find(endsWith(fnames,'_RasterPlot.fig'),14)) && ~isempty(find(strcmp(fnames,'ms.mat'),6)) && ~isempty(find(strcmp(fnames,'SFP.mat'),6))% && ~isempty(find(contains(folders{i},'Miniscope'),1));
            cd(folders{i});    
            a = pwd;
            msfile = append(pwd,'\ms.mat');
            splitends = strsplit(a,'\');
            
            ID = splitends(end); 
            
            try
                load(msfile)
                NumNeur = num2str(ms.numNeurons);
                if show == false
                    f = figure('visible','off');
                    set(gca,'FontSize',3);
                else
                    f = figure('visible','on');
                    set(gca,'FontSize',6);
                end
                if save == true
                    set(gca,'FontSize',3);
                end
                subplot(2,1,1)
                imagesc(max(ms.SFPs,[],3))
                title([ID '--' 'N = ' NumNeur], 'Interpreter', 'none')
                pbaspect([1 1 1])
                subplot(2,1,2)
                imagesc(ms.FiltTraces')
                title('Calcium Trace Raster')
                pbaspect([1 1 1])
                colorbar
                caxis([cmin cmax])
                if save == true
                    savefig([ID '_RasterPlot.fig'])
                    saveas(f,[ID '_RasterPlot' fileformat])
                    display(['saved ' ID])
                    pause(1)
                end
            catch
                display(['Failed to generate figure for ',ID])
                continue
            end            
        end
    end
end
end