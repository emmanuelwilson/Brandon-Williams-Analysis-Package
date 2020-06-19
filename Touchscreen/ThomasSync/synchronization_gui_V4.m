function synchronization_gui_V4()

% Function that reads the timestamp.csv, schedules.csv files and videos
% labeled msCam[].avi and behaviorCam[].avi, and builds synchronization
% matrices that are saved in msTouchSynch.mat.

% timestamp.csv: file produced by the experiment containing times when
% cameras 0 and 1 (measured by the system clock) took frames, which were
% then saved in video files.

% schedules.csv: file that contains the output from the cage during the
% recording session, with time measured relative to the cage clock.

% ---------------- display GUI to select dataset --------------------------

% read screen size
screensize = get(0,'Screensize');

% flag to skip datasets if there are problems
dataset_ok = 0;

% parameters of the task
complexity = 'STAGE 1';
separation = 'S3';
delay = '6 SEC';

% list of subfolders to verify
list_subfolders = '';

% 0 = folders under tunl/stage 1/etc., 1 = specific folder, 2 = all subfolders
% below a given folder
subfolder = 0;

% build the GUI here
gui_fig = figure('Position',[screensize(3)/3 screensize(4)/2 screensize(3)/3 screensize(4)/6],...
    'Name','Select the folder containing your data...','NumberTitle','off');

my_data_dir = pwd;

% text box
txtbox = uicontrol('Style','edit',...
    'Position',[screensize(4)/40 screensize(3)/20 screensize(4)/3 screensize(3)/50],...
    'String',my_data_dir);
txtbox.Units = 'normalized';

% browse button
browse_but = uicontrol('Style','pushbutton',...
    'Position',[0.4*screensize(4) screensize(3)/20 0.1*screensize(4) screensize(3)/50],...
    'String','Browse','Callback',@f_browse);
browse_but.Units = 'normalized';

% task complexity (stage 1, 2 and 3)
txt = uicontrol('Style','text',...
    'Position',[screensize(4)*0.0125 0.02*screensize(3) 0.075*screensize(4) screensize(3)/50],...
    'String','task complexity');
txt.Units = 'normalized';

txtcomplexity = uicontrol('Style','edit',...
    'Position',[screensize(4)/40 0.01*screensize(3) 0.05*screensize(4) screensize(3)/50],...
    'String',complexity);
txtcomplexity.Units = 'normalized';

% stimulus separation (S1, 2 and 3)
txt1 = uicontrol('Style','text',...
    'Position',[screensize(4)*(0.015+0.1) 0.02*screensize(3) 0.075*screensize(4) screensize(3)/50],...
    'String','stimulus separation');
txt1.Units = 'normalized';

txtseparation = uicontrol('Style','edit',...
    'Position',[screensize(4)*(1/40+0.1) 0.01*screensize(3) 0.05*screensize(4) screensize(3)/50],...
    'String',separation);
txtseparation.Units = 'normalized';

% delay (2, 4 and 6 SEC)
txt3 = uicontrol('Style','text',...
    'Position',[screensize(4)*(1/40+2*0.1) 0.02*screensize(3) 0.05*screensize(4) screensize(3)/50],...
    'String','delay duration');
txt3.Units = 'normalized';

txtdelay = uicontrol('Style','edit',...
    'Position',[screensize(4)*(1/40+2*0.1) 0.01*screensize(3) 0.05*screensize(4) screensize(3)/50],...
    'String',delay);
txtdelay.Units = 'normalized';

% checkbox to select if we want to check a folder or a subfolder
checkfolder = uicontrol('Style','checkbox',...
    'Position',[screensize(4)*(1/40+3*0.1) 0.01*screensize(3) 0.05*screensize(4) screensize(3)/50],...
    'String','subfolder','Value',0,'Callback',@f_check);
checkfolder.Units = 'normalized';

checkfolders = uicontrol('Style','checkbox',...
    'Position',[screensize(4)*(1/40+3*0.1) 0.03*screensize(3) 0.05*screensize(4) screensize(3)/50],...
    'String','subfolders','Value',0,'Callback',@f_checks);
checkfolders.Units = 'normalized';

% proceed button
proceed_but = uicontrol('Style','pushbutton',...
    'Position',[0.4*screensize(4) 0.01*screensize(3) 0.1*screensize(4) screensize(3)/50],...
    'String','Proceed','Callback',@f_proceed);
proceed_but.Units = 'normalized';

% ================================= functions =============================

% gui functions

% checkbox specifying if the user wants to analyse a single recording
% session, or a bunch contained in a folder
    function f_check(source,event)
        
        if source.Value==1
            % analysis of a single folder. Erase the complexity, separation
            % and delay
            txtcomplexity.String = '';
            txtseparation.String = '';
            txtdelay.String = '';
            % just the specified folder
            subfolder = 1;
            
            % reset subfolders box
            set(checkfolders,'Value',0);
        else
            % analysis of a folder containing several recording sessions
            txtcomplexity.String = complexity;
            txtseparation.String = separation;
            txtdelay.String = delay;
            % all subfolders specified by complexity, etc.
            subfolder = 0;
        end
        
    end

% checkbox specifying if the user wants to analyse all folders immediately
% below a given folder
    function f_checks(source,event)
        
        if source.Value==1
            % analysis of a single folder. Erase the complexity, separation
            % and delay
            txtcomplexity.String = '';
            txtseparation.String = '';
            txtdelay.String = '';
            % all subfolders below a given folder
            subfolder = 2;
            
            % reset subfolder box
            set(checkfolder,'Value',0);
        else
            % analysis of a folder containing several recording sessions
            txtcomplexity.String = complexity;
            txtseparation.String = separation;
            txtdelay.String = delay;
            % all subfolders specified by complexity, etc.
            subfolder = 0;
        end
        
    end

% point to the folder containing the data
    function f_browse(source,event)
        
        temp = uigetdir(my_data_dir,'Please select folder');
        if temp~=0
            my_data_dir = temp;
        end
        set(txtbox,'String',my_data_dir);
        
        % show content of selected folder
        switch subfolder
            
            case 0
                % we are checking a whole set of sessions (tunl/stage 1/etc.) contained in folder
                my_dir = [my_data_dir '\TUNL\' complexity '\' separation '\' delay];
                temp = dir(my_dir);
                
                list_subfolders = [];
                for i1=1:length(temp)
                    if strcmp(temp(i1).name,'.')==0 && ...
                            strcmp(temp(i1).name,'..')==0
                        list_subfolders{end+1} = [my_data_dir '\TUNL\' complexity '\' separation '\' delay '\' temp(i1).name];
                    end
                end
                
                % display what has been found
                msgbox(list_subfolders)
                
            case 1
                % we are checking a single folder which we point to precisely
                my_dir = my_data_dir;
                
                % display what has been found
                msgbox(my_dir)
            case 2
                % we are checking a whole set of subfolders contained in folder
                my_dir = my_data_dir;
                temp = dir(my_dir);
                
                list_subfolders = [];
                for i1=1:length(temp)
                    % check that is is a folder
                    if strcmp(temp(i1).name,'.')==0 && ...
                            strcmp(temp(i1).name,'..')==0 && ...
                            temp(i1).isdir
                        list_subfolders{end+1} = [my_data_dir '\' temp(i1).name];
                    end
                end
                
                % display what has been found
                msgbox(list_subfolders)
        end
        
    end

    function f_proceed(source,event)
        
        disp('Starting processing data...');
        close(gui_fig);
        
        % User can either select directly the subfolder containing the
        % videos, etc., or a higher folder containing several such
        % subfolders
        
        if subfolder==1
            % we are checking a single folder which we point to precisely
            my_path = my_data_dir;
            dataset_ok = 0;
            
            disp(' ');
            disp('================================================================================================');
            disp(['Folder ' my_path]);
            str = '';
            for j1=1:length(my_path)
                str = [str '--------'];
            end
            disp(str);
            disp(' ');
            
            % check if timestamp.csv file exists
            resp = exist([my_path '\timestamp.csv'],'file');
            if resp==0
                waitfor(msgbox(['Cannot find file ' my_path '\timestamp.csv. Stopping here...']));
            else
                dataset_ok = dataset_ok + 1;
            end
            
            % check if schedules.csv file exists
            resp = exist([my_path '\schedules.csv'],'file');
            if resp==0
                waitfor(msgbox(['Cannot find file ' my_path '\schedules.csv. Stopping here...']));
            else
                dataset_ok = dataset_ok + 1;
            end
            
            % check if msCam[].avi files exist
            resp = exist([my_path '\msCam1.avi'],'file');
            if resp==0
                waitfor(msgbox(['Cannot find file ' my_path '\msCam1.avi. Stopping here...']));
            else
                dataset_ok = dataset_ok + 1;
            end
            
            % check if behavCam[].avi files exit
            resp = exist([my_path '\behavCam1.avi'],'file');
            if resp==0
                waitfor(msgbox(['Cannot find file ' my_path '\behavCam1.avi. Stopping here...']));
            else
                dataset_ok = dataset_ok + 1;
            end
            
            if dataset_ok==4
                single_folder_synchronization(my_path);
                close all
            end
        else
            % we are checking a whole set of sessions contained in folder
            % (either tunl/stage 1/etc. or subfolders in a given folder).
            
            % process one-by-one all folders found in the main directory
            for i1=1:length(list_subfolders)
                
                dataset_ok = 0;
                
                my_path = char(list_subfolders(i1));
                
                disp(' ');
                disp('================================================================================================');
                disp(['Folder ' my_path]);
                str = '';
                for j1=1:length(my_path)
                    str = [str '--------'];
                end
                disp(str);
                disp(' ');
                
                % check if timestamp.csv file exists
                resp = exist([my_path '\timestamp.csv'],'file');
                if resp==0
                    waitfor(msgbox(['Cannot find file ' my_path '\timestamp.csv. Stopping here...']));
                else
                    dataset_ok = dataset_ok + 1;
                end
                
                % check if schedules.csv file exists
                resp = exist([my_path '\schedules.csv'],'file');
                if resp==0
                    waitfor(msgbox(['Cannot find file ' my_path '\schedules.csv. Stopping here...']));
                else
                    dataset_ok = dataset_ok + 1;
                end
                
                % check if msCam[].avi files exist
                resp = exist([my_path '\msCam1.avi'],'file');
                if resp==0
                    waitfor(msgbox(['Cannot find file ' my_path '\msCam1.avi. Stopping here...']));
                else
                    dataset_ok = dataset_ok + 1;
                end
                
                % check if behavCam[].avi files exit
                resp = exist([my_path '\behavCam1.avi'],'file');
                if resp==0
                    waitfor(msgbox(['Cannot find file ' my_path '\behavCam1.avi. Stopping here...']));
                else
                    dataset_ok = dataset_ok + 1;
                end
                
                if dataset_ok==4
                    single_folder_synchronization(my_path);
                    close all
                end
                
                disp(' ');
                disp('================================================================================================');
                disp(' ');
                
            end
        end
        
    end

end






















