%Author : Tori-Lynn Temple, A part of Spatial Coding for JB and Guillaume 

%This function is coding for the binning array and outputs per frame the
%appropriate bin location for each individual cell.
%This is done by taking the location from ms and behav of each cell
%and converting to the location in cm of each cell. From here we set a
%layout for the bins, meaning each bin is at a certain location in the
%field space in cm, and we place the cells into each. The out put of this
%function is then the binning array at a single frame. 

function [ out ] = BinningTheCells(ms, behav, first_frame)


%  TIME AND FRAME CORRECTION IN ORDER TO FIND CELL LOCATIONS
%This is to make sure every time vector has only unique values
 
[unique_behavtime, IAbehav, ICbehav]=unique(behav.time);
[unique_mstime, IAms, ICms]=unique(ms.time);

% Extract one trace of a give cell_ID cell, excluding doubled timestamps:
trace = ms.trace(IAms,cell_ID);

% Interpolate the X position of the mouse. 
interpolated_X=interp1(behav.time(IAbehav),behav.position(IAbehav,1),ms.time(IAms),'pchip');
interpolated_Y=interp1(behav.time(IAbehav),behav.position(IAbehav,2),ms.time(IAms),'pchip');

%%  USER INPUT VARIABLES  

%The use inputs the sizw of the field so the program can properly form a
%well sized mesh of bins in the field and assign cells appropriatley. 
prompt = 'Enter in the width of the field in centimeters:';
length_of_field_x = input(prompt);

prompt = 'Enter in the length of the field in centimeters:';
length_of_field_y = input(prompt);

prompt = 'Enter in the size of the bins you prefer, NORMALLY PEOPLE CHOOSE 5CM:';
bin_sizes = input(prompt);

prompt = 'Finally click and drag the mouse to measure the x-field in pixels:';
figure, imshow(first_frame);  
h = imline; 
drawn_line = api.getPosition(h); %where draw_line becomes a 2x2 array of [X1 Y1; X2 Y2]
draw_line_x_length = draw_line(2,1) -draw_line(1,1) ; %returns the length of line in pixels 

%%  CREATING THE BIN ARRAY BASED ON VARIABLE INPUT

%now that we can find the area of the field, we apply an appropriate
%number of bins. 
%where bins_x = the number of bins in the x directions 
%and bins_y = the number of bins in the y direcitons 
bins_x = length_of_field_x/bin_sizes;  
bins_y = length_of_field_y/bin_sizes;  

%setting the bin array based on the size of the field
cell_bins = zeros(bins_y,bins_x) ; 


%% FILLING THE BIN ARRAY WITH CELLS IN EACH FRAME 
%now that the bin array has been made we can fill the cells into the
%appropriate bins by a simple algorith below (explanaition of the 
%algorithm can be found in the document written about this script on 
%github. See:________________________):

%first setting a vector space for all the cells positions. 
i = 1:1:length(ms.cells);

%a pixel to cm conversion for the drawn line. This makes it so we can
%compare this with the entered field in cm with the drawn field in pixels.
%Which we will use a conversion for later calculations. 
pixels_to_cm = length_of_field_x/drawn_line ; 

%converts the found positions of x and y of a cell to cm
location_of_cell_X(i) = interpolated_X(i)*pixels_to_cm; %In Centimeters 
location_of_cell_Y(i) = interpolated_Y(i)*pixels_to_cm; %In Centimeters 

 %Now we use the x and y positions in cm to calculate the right bin placement.  
 bin_cell_location_x(i) = ceil(location_of_cell_X(i)/bin_sizes); 
 bin_cell_location_y(i) = ceil(location_of_cell_Y(i)/bin_sizes);  
 
 
 %once it finds the proper bin to place the cell, we place it in the bin
 %array as the proper location for the cell. 
 cell_bins = (bin_cell_location_y(i) , bin_cell_location_x(i)); 

end 