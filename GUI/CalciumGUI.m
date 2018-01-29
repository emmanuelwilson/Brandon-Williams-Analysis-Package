function varargout = CalciumGUI(varargin)
% CALCIUMGUI MATLAB code for CalciumGUI.fig
%      CALCIUMGUI, by itself, creates a new CALCIUMGUI or raises the existing
%      singleton*.
%
%      H = CALCIUMGUI returns the handle to a new CALCIUMGUI or the handle to
%      the existing singleton*.
%
%      CALCIUMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALCIUMGUI.M with the given input arguments.
%
%      CALCIUMGUI('Property','Value',...) creates a new CALCIUMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CalciumGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CalciumGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CalciumGUI

% Last Modified by GUIDE v2.5 16-Jan-2018 14:15:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CalciumGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CalciumGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This GUI is intended to help with exploring and providing a greater
%visual understanding of neural data. This particular itteration provides a
%cell by cell look of the calcium fluorescence trace, a trajectory plot,
%polar plot and location where it fired as well as a physiological look of 
%the cell in question. This GUI will also provide with a firing option, 
%which looks at the estimated cell firing based on the fluorescence. 
%PLEASE PROVIDE FILE NAMES IN LINES 69-71 TO BE PROCESSED/LOOKED AT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Author Emmanuel Wilson

% --- Executes just before CalciumGUI is made visible.
function CalciumGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CalciumGUI (see VARARGIN)

% Choose default command line output for CalciumGUI
handles.output = hObject;
handles.cell = 1;                   %Cell initialization
handles.thresh = 0.1;               %Threshold default

[trackName, trackPath, trackFilter] = uigetfile('*.mat','Trajectory Head Tracking Data');
[HeadName, HeadPath, HeadFilter] = uigetfile('*.mat','Head Direction Data');
[frameName, framePath, frameFilter] = uigetfile('*.mat','Frame Map');
[msName, msPath, msFilter] = uigetfile('*.mat','Physiology struct .mat file');
trackLoc = strcat(trackPath, trackName);
HeadLoc = strcat(HeadPath, HeadName);
frameLoc = strcat(framePath, frameName);
msLoc = strcat(msPath, msName);
a = load(trackLoc);
b = load(HeadLoc);
c= load(frameLoc);
d= load(msLoc);
handles.tracking = a.SINKdata;      %Accessing Head location trajectory map
handles.HD = b.HDdeg;               %Accessing Head direction 
handles.sinkedframes = c.frameMap;  %Accessing frame map
handles.ms = d.ms;                  %Accessing ms
handles.ms.firing_cnmfe = handles.ms.firing_cnmfe; %modify cnmfe firing for visual ease
handles.riseT = 10;                 %Default rise time
handles.frame = length(handles.sinkedframes);   %# of frames in videos
handles.lframe = 1;                 %Lower frame bounds (default 0)
handles.firingON = 0;               %Active if firing trajectory map is on
handles.firingTraceVogle = 0;       %Active if Voglestein trace is on
handles.firingTraceCnmfe = 0;        %Active if Pnevmatikakis trace is on
handles.htrack = [];
guidata(hObject, handles);          %save parameter changes
updateAxes(handles)                 %update graphs


% --- Outputs from this function are returned to the command line.
function varargout = CalciumGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function ThreshSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ThreshSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject, 'Min',0);                                  %set minimum value                        
set(hObject, 'Max',0.25) ;                              %set maximum value
uicontrol('Style','slider','Min',0,'Max',0.25,'Value',0.1)   %UI slider interface
handles.thresh = get(hObject, 'Value');                 %Set threshold to slider value
guidata(hObject,handles)                                %save changes
updateAxes(handles)                                     %Update axes


% --- Executes during object creation, after setting all properties.
function ThreshSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThreshSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function TimeSlide_Callback(hObject, eventdata, handles)
% hObject    handle to TimeSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject, 'Min',0);                                                          %set minimum value
set(hObject, 'Max',length(handles.sinkedframes)) ;                              %set maximum value    
uicontrol('Style','slider','Min',1,'Max',length(handles.sinkedframes),'Value',1)%UI slider interface
handles.frame = get(hObject, 'Value');                                          %set frame number to slider value    
guidata(hObject,handles)                                                        %Save changes 
updateAxes(handles)                                                             %Update axes


% --- Executes during object creation, after setting all properties.
function TimeSlide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function TimeSlideDown_Callback(hObject, eventdata, handles)
% hObject    handle to TimeSlideDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject, 'Min',0);                                                          %set minimum value
set(hObject, 'Max',length(handles.sinkedframes)) ;                              %set maximum value
uicontrol('Style','slider','Min',0,'Max',length(handles.sinkedframes),'Value',10)%UI slider interface
if get(hObject, 'Value')>1
    handles.lframe = get(hObject, 'Value');                                         %set frame number to slider value
else
    handles.lframe = 1;
end
guidata(hObject,handles)                                                        %Save changes 
updateAxes(handles)                                                             %update axes


% --- Executes during object creation, after setting all properties.
function TimeSlideDown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeSlideDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in BackButton.
function BackButton_Callback(hObject, eventdata, handles)
% hObject    handle to BackButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.cell) || handles.cell == 1       %Back Button actions 
    handles.cell = 1;                               %default to cell#1 & set limit to avoid going below cell#1
else
    handles.cell = handles.cell - 1;                %Go back one cell
end
guidata(hObject,handles)                            %save changes
updateAxes(handles)                                 %update axes


% --- Executes on button press in Forwardbutton.
function Forwardbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Forwardbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%FORWARD button actions
if isempty(handles.cell)
    handles.cell = 1;                                   %default cell#1
elseif handles.cell >= length(handles.ms.trace(1,:))    
    handles.cell = length(handles.ms.trace(1,:));       %set max cell# within the data provided
else
    handles.cell = handles.cell + 1;                    %Go forward one cell
end
guidata(hObject,handles)                %save changes 
updateAxes(handles)                     %update axes

function updateAxes(handles)
set(handles.CellNumView, 'String', num2str(handles.cell))                  %Set Cell # display window 
set(handles.ThreshView, 'String', num2str(handles.thresh))                 %Set threshold display
set(handles.FrameNumView, 'String', num2str(handles.frame))                %Set upper bound frame # display
set(handles.FrameNumViewDown, 'String',num2str(handles.lframe))            %Set lower bound frame # display
trace = handles.ms.trace(:,handles.cell);                                  %Extract trace
if handles.sinkedframes(length(handles.sinkedframes)) < length(handles.ms.trace)    
    trace = trace(handles.sinkedframes);                                   %Sink trace if not already sinked
end
  
shifted_trace = trace((handles.riseT : length(trace)));                    %Shift trace by time offset riseT
shifted_position_track = (handles.tracking(1 : (length(handles.tracking) - handles.riseT + 1)))';   %track position with offset
shifted_trace(shifted_trace < handles.thresh) = 0;                         %apply threshold on tracking
ind_firing_on_trajectory = find(shifted_trace);                            %find all relevant points for head trajectory

if handles.firingTraceVogle == 1
    firing = handles.ms.firing(:,handles.cell);                            %Extract firing trace
    fire = firing;                                                         %Duplicate
    fire(fire < handles.thresh) = 0;                                       %apply threshold on firing
    for i = 1 : length(trace)                                              %Apply frame bounds on graphs
        if i > handles.frame || i < handles.lframe
            shifted_trace(i) = 0;
            fire(i) = 0;
        end
    end
    ind_fire = find(fire);                                                 %find all relevant points for firing trajectory
end
if handles.firingTraceCnmfe == 1
    firingc = handles.ms.firing_cnmfe(:,handles.cell);                     %Extract firing trace
    firec = firingc;                                                       %Duplicate
    firec(firec < handles.thresh) = 0;                                     %apply threshold on firing
    for i = 1 : length(trace)                                              %Apply frame bounds on graphs
        if i > handles.frame || i < handles.lframe
            shifted_trace(i) = 0;
            firec(i) = 0;
        end
    end
    ind_firec = find(firec);                                               %find all relevant points for firing trajectory   
end

if handles.firingON == 1 && handles.firingTraceVogle == 1 && handles.firingTraceCnmfe == 0
   htrack(:,1) = handles.tracking(handles.sinkedframes(ind_fire),1);       %change coordinate system when Voglestein firing is turned ON and CNMMFE firing is OFF
   htrack(:,2) = handles.tracking(handles.sinkedframes(ind_fire),2);
   headdir = handles.HD(handles.sinkedframes(ind_fire));                   %set head direction matrix to correspond with current firing mode
elseif handles.firingON == 1 && handles.firingTraceCnmfe == 1 && handles.firingTraceVogle == 0
    htrack(:,1) = handles.tracking(handles.sinkedframes(ind_firec),1);     %change coordinate system when CNMFE Firing is turned ON AND Voglestein is OFF
    htrack(:,2) = handles.tracking(handles.sinkedframes(ind_firec),2);
    headdir = handles.HD(handles.sinkedframes(ind_firec));                 %set head direction matrix to correspond with current firing mode
else
    htrack(:,1) = handles.tracking(handles.sinkedframes(ind_firing_on_trajectory),1);   %Set coordinate system when firing is OFF or multiple firing modes are ON
    htrack(:,2) = handles.tracking(handles.sinkedframes(ind_firing_on_trajectory),2);
    headdir = handles.HD(handles.sinkedframes(ind_firing_on_trajectory));  %set head direction matrix to fluorescence
end
%calcium fluorescence plot
axes(handles.axes1)
if handles.firingTraceVogle == 0 && handles.firingTraceCnmfe == 0
    plot(trace,'b')
elseif handles.firingTraceVogle == 1 && handles.firingTraceCnmfe == 1
    plot(trace,'b')
    hold on
    plot(firing,'m')
    plot(firingc,'r')
    hold off
elseif handles.firingTraceVogle == 1
    plot(trace,'b')
    hold on
    plot(firing,'m')
    hold off
elseif handles.firingTraceCnmfe == 1
    plot(trace,'b')
    hold on
    plot(firingc,'r')
    hold off
end

hline(handles.thresh)
vline(handles.frame)
vline(handles.lframe)
title('Calcium Transiance')
xlabel('Time(frames)')
ylabel('Amplitude')

%Trajectory plot
axes(handles.axes2)
plot(handles.tracking(handles.sinkedframes(int16(handles.lframe):int16(handles.frame)),1),-handles.tracking(handles.sinkedframes(int16(handles.lframe):int16(handles.frame)),2),'k')
 
colormap(hsv)
title('Trajectory Map')
xlabel('Position in X(Pixel)')
ylabel('Position in Y(Pixel)')
hold on
scatter(htrack(:,1),-htrack(:,2),50,headdir,'filled')
hold off
caxis([0 360])

%Head direction trajectory legend
axes(handles.axes4)
r = 0.5:0.5:1;
theta = linspace(0, 2*pi, 100);
[rg, thg] = meshgrid(r,theta);
[x,y] = pol2cart(thg,rg);
p = pcolor(x,-y,thg);
camroll(90)
colormap(hsv);
shading flat;
axis off;
%physiology plot
axes(handles.axes3)
imshow(uint8(handles.ms.frameMax) + 255*uint8(bwperim(handles.ms.segments(:,:,handles.cell)))) % replace ms.meanFrame with ms.frameMax
title('Physiology- Cell position')
%polar plot
[sorted_HD_no_shift, ind_sorted_HD_no_shift] = sort(deg2rad(headdir));
fire_temp = handles.ms.firing(ind_sorted_HD_no_shift);
fire_temp = fire_temp';
axes(handles.axes5)
binned_sorted_HD_no_shift = [];
binned_firing = [];
sorted_HD_no_shift = floor(rad2deg(sorted_HD_no_shift));
sorted_HD_no_shift = sorted_HD_no_shift';
ang_bin_size = 5;
i = 1;
temp_sorted_HD_no_shift = [];
temp_firing = [];
for ind = 1 : 360/ang_bin_size
    k = 1;
    while (i < (length(sorted_HD_no_shift) + 1)) && (sorted_HD_no_shift(i) < (ang_bin_size*ind + sorted_HD_no_shift(1)))
        temp_sorted_HD_no_shift(k) = sorted_HD_no_shift(i);
        temp_firing(k) = fire_temp(i);% handles.ms.firing(i,handles.cell);
        i = i + 1;
        k = k + 1;
    end
    binned_sorted_HD_no_shift(ind) = mean(temp_sorted_HD_no_shift);
    binned_firing(ind) = mean(temp_firing);
    temp_sorted_HD_no_shift = [];
    temp_firing = [];
end

binned_sorted_HD_no_shift = deg2rad(binned_sorted_HD_no_shift);
binned_sorted_HD_no_shift(isnan(binned_sorted_HD_no_shift)) = [];
binned_firing(isnan(binned_firing)) = [];
polar(binned_sorted_HD_no_shift, ((binned_firing)));

title('Polar Plot of Head Direction During Active Period')


% --- Executes on button press in Vogle.
function Vogle_Callback(hObject, eventdata, handles)
% hObject    handle to Vogle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.firingTraceVogle = get(hObject,'Value'); %boxes checked yes or no
guidata(hObject,handles)                    %save changes 
updateAxes(handles)                         %update axes


% --- Executes on button press in firingON.
function firingON_Callback(hObject, eventdata, handles)
% hObject    handle to firingON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.firingON = get(hObject,'Value');    %Box checked, yes  or no
guidata(hObject,handles)                    %save changes
updateAxes(handles)                         %update axes



function CellNumView_Callback(hObject, eventdata, handles)
% hObject    handle to CellNumView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

celVar = get(handles.CellNumView, 'string');            %Get string from text box
celVar = str2double(celVar);                            %convert to double
if celVar > 0 && celVar <= length(handles.ms.trace(1,:))     %update values
    handles.cell = celVar;    
end
guidata(hObject,handles)                                %Save changes
updateAxes(handles)                                     %update axes


% --- Executes during object creation, after setting all properties.
function CellNumView_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellNumView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ThreshView_Callback(hObject, eventdata, handles)
% hObject    handle to ThreshView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

threshVar = get(handles.ThreshView, 'string');      %Get string from text box
threshVar = str2double(threshVar);                  %convert to double
if threshVar > 0 && threshVar <= 0.25               %update values
    handles.thresh = threshVar;
end
guidata(hObject,handles)                            %Save changes
updateAxes(handles)                                 %update axes


% --- Executes during object creation, after setting all properties.
function ThreshView_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThreshView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FrameNumViewDown_Callback(hObject, eventdata, handles)
% hObject    handle to FrameNumViewDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DownFrameVar = get(handles.FrameNumViewDown, 'string');               %Get string value from text box
DownFrameVar = str2double(DownFrameVar);                            %convert to double
if DownFrameVar > 0 && DownFrameVar <= length(handles.tracking)     %Update values
    handles.lframe = DownFrameVar;
end
guidata(hObject,handles)                                        %Save changes
updateAxes(handles)                                             %Update axes



% --- Executes during object creation, after setting all properties.
function FrameNumViewDown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameNumViewDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FrameNumView_Callback(hObject, eventdata, handles)
% hObject    handle to FrameNumView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
upFrameVar = get(handles.FrameNumView, 'string');               %Get string value from text box
upFrameVar = str2double(upFrameVar);                            %convert to double
if upFrameVar > 0 && upFrameVar <= length(handles.tracking)     %Update values
    handles.frame = upFrameVar;
end
guidata(hObject,handles)                                        %Save changes
updateAxes(handles)                                             %Update axes


% --- Executes during object creation, after setting all properties.
function FrameNumView_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameNumView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CNMFEcheck.
function CNMFEcheck_Callback(hObject, eventdata, handles)
% hObject    handle to CNMFEcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.firingTraceCnmfe = get(hObject,'Value'); %boxes checked yes or no
guidata(hObject,handles)                    %save changes 
updateAxes(handles)                         %update axes
% Hint: get(hObject,'Value') returns toggle state of CNMFEcheck
