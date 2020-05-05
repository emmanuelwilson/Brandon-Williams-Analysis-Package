function [ Cell location in time  ] = cell_location(ms, behav, cell, time)

%% Parameters
min_speed = 2;
smoothing = 1;

%% Extracting data from structure
calcium_time = ms.time;
dt = mode(diff(calcium_time));  %returns the frequent element in the calcium_time array after taking the difference between each element 
Fs = 1/dt;

%% 

%%Changing the forloop to indexing: 
% i = 1:1:length(ms.time); 
% behav_time_at_calcium_time_idx(i) = dsearchn(behav.time,ms.time(i));

%Re-aligning calcium and behavior times
for i = 1:length(ms.time)
behav_time_at_calcium_time_idx(i) = dsearchn(behav.time,ms.time(i));
end

X_at_calcium_time = behav.position(behav_time_at_calcium_time_idx,1);
Y_at_calcium_time = behav.position(behav_time_at_calcium_time_idx,2);
speed_at_calcium_time = behav.speed(behav_time_at_calcium_time_idx);


position_at_calcium_time(:,1) = X_at_calcium_time;
position_at_calcium_time(:,2) = Y_at_calcium_time;

low_speed_idx = find(speed_at_calcium_time<min_speed);



end 