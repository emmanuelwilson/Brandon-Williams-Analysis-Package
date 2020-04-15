function ms = mkStablizedVideos(doPath)

%     cleanCalciumVideos(doPath)

    %% msRun2018
    % Version 1.0 GE
    % Updated version of the msRun script originally proposed by Daniel B
    % Aharoni to analyse miniscope 1p calcium imaging data.
    % This version is build on top of the original package to maximize compatibility.
    % It includes NormCorre for image registration, CNMF-E for source extraction,
    % and CellReg for chronic registration across sessions. It also includes
    % custom written scripts to explore the data (eg spatial firing, transients
    % properties visualization)

    % Copyright (C) 2017-2018 by Guillaume Etter
    %
    % This program is free software; you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation; either version 3 of the License, or any
    % later version.
    % Contact: etterguillaume@gmail.com

    %% Auto-detect operating system
    if ispc
        separator = '\'; % For pc operating systems
    else
        separator = '/'; % For unix (mac, linux) operating systems
    end
    
    
    %% Parameters
    spatial_downsampling = 3; % (Recommended range: 2 - 4. Downsampling significantly increases computational speed, but verify it does not
    isnonrigid = false; % If true, performs non-rigid registration (slower). If false, rigid alignment (faster).
    analyse_behavior = false;
    copy_to_googledrive = false;
    if copy_to_googledrive;
        copydirpath = uigetdir([],'Please select the root folder in which files will be copied');
    end

    % Generate timestamp to save analysis
    script_start = tic;
    analysis_time =strcat(date,'_', num2str(hour(now)),'-',num2str(minute(now)),'-',num2str(floor(second(now))));

    %% 1 - Create video object and save into matfile
    fprintf('\n\tStep 1: Create video object');
    ms = msGenerateVideoObj(doPath,'msCam');
    ms.analysis_time = analysis_time;
    ms.ds = spatial_downsampling;
    ms.dirName = strcat(doPath,separator,analysis_time);
    mkdir(ms.dirName);
    
    %% 2 - Perform motion correction using NormCorre
    
    fprintf('\n\tStep 2: Motion correction\n');
    ms = msNormCorre(ms,isnonrigid,['StablizedVideos/' doPath(find(ismember(doPath,'/'),1,'first')+1:end)]);
end