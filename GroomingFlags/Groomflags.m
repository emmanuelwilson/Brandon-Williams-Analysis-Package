function varargout = Groomflags(varargin)
% GROOMFLAGS MATLAB code for Groomflags.fig
%      GROOMFLAGS, by itself, creates a new GROOMFLAGS or raises the existing
%      singleton*.
%
%      H = GROOMFLAGS returns the handle to a new GROOMFLAGS or the handle to
%      the existing singleton*.
%
%      GROOMFLAGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GROOMFLAGS.M with the given input arguments.
%
%      GROOMFLAGS('Property','Value',...) creates a new GROOMFLAGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Groomflags_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Groomflags_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Groomflags

% Last Modified by GUIDE v2.5 06-Sep-2018 14:27:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Groomflags_OpeningFcn, ...
                   'gui_OutputFcn',  @Groomflags_OutputFcn, ...
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


% --- Executes just before Groomflags is made visible.
function Groomflags_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Groomflags (see VARARGIN)

% Choose default command line output for Groomflags
handles.output = hObject;
handles.fpath = uigetdir;
handles.vidframe = 1;
handles.vidname = '';
handles.initialized = 0;
handles.grooming = zeros(1,2);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Groomflags wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Groomflags_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function FrameSlide_Callback(hObject, eventdata, handles)
% hObject    handle to FrameSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject, 'Min',1);                                  %set minimum value
set(hObject, 'Max',1000) ;                              %set maximum value    
uicontrol('Style','slider','Min',1,'Max',1000,'Value',1)%UI slider interface
handles.vidframe = round(get(hObject, 'Value'));                                          %set frame number to slider value    
guidata(hObject,handles)                           %Save changes 
updateAxes(handles)     


% --- Executes during object creation, after setting all properties.
function FrameSlide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function FrameDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to FrameDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FrameVar = get(handles.FrameDisplay, 'string');               %Get string value from text box
FrameVar = str2double(FrameVar);                            %convert to double
if FrameVar > 0 && FrameVar <= handles.vid.Duration*handles.vid.FrameRate    %Update values
    handles.vidframe = FrameVar;
end
guidata(hObject,handles)                                        %Save changes
updateAxes(handles)                                             %Update axes


% --- Executes during object creation, after setting all properties.
function FrameDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Vidlist.
function Vidlist_Callback(hObject, eventdata, handles)
% hObject    handle to Vidlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.contents = dir([handles.fpath '\behav*.avi']);
vlist = {handles.contents(~[handles.contents.isdir]).name};
set(handles.Vidlist,'String',vlist)
index = get(handles.Vidlist,'Value');
filelist = get(handles.Vidlist,'String');
handles.selected = filelist{index};
handles.vidname = [handles.fpath '\' filelist{index}];
handles.vid = VideoReader(handles.vidname);

guidata(hObject,handles)
updateAxes(handles)

% Hints: contents = cellstr(get(hObject,'String')) returns Vidlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Vidlist


% --- Executes during object creation, after setting all properties.
function Vidlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vidlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in playpausebutton.
function playpausebutton_Callback(hObject, eventdata, handles)
% hObject    handle to playpausebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% vid = VideoReader(handles.vidname);
axes(handles.axes1)
while get(hObject,'Value')   
    image(read(handles.vid,handles.vidframe))
    if handles.vidframe < handles.vid.Duration*handles.vid.FrameRate
        handles.vidframe = handles.vidframe +1;    
        guidata(hObject,handles)
    else
        handles.vidframe = 1;
        guidata(hObject,handles)
        break
    end
    pause(0.05)
end
guidata(hObject,handles)
updateAxes(handles)


% --- Executes on button press in FlagStop.
function FlagStop_Callback(hObject, eventdata, handles)
% hObject    handle to FlagStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.grooming(end,2) == 0
    handles.grooming(end,2) = handles.vidframe + (str2num(handles.selected(9:end-4))-1)*1000;
else
    handles.grooming(end+1,2) = handles.vidframe + (str2num(handles.selected(9:end-4))-1)*1000;
end
guidata(hObject,handles)
updateAxes(handles)

% --- Executes on button press in FlagStart.
function FlagStart_Callback(hObject, eventdata, handles)
% hObject    handle to FlagStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.grooming(1,1) == 0
    handles.grooming(1,1) = handles.vidframe + (str2num(handles.selected(9:end-4))-1)*1000;
else
    handles.grooming(end+1,1) = handles.vidframe + (str2num(handles.selected(9:end-4))-1)*1000;
end

guidata(hObject,handles)
updateAxes(handles)

function updateAxes(handles)
% if handles.newvid
%     handles.vid = VideoReader(handles.vidname);
%     handles.newvid = 0;
% end
x = (str2num(handles.selected(9:end-4))-1)*1000:1:(str2num(handles.selected(9:end-4))-1)*1000+ 1000;
xmin = (str2num(handles.selected(9:end-4))-1)*1000;
xmax = (str2num(handles.selected(9:end-4))-1)*1000+ 1000;
set(handles.FrameDisplay, 'String', num2str(handles.vidframe))
axes(handles.axes1)
image(read(handles.vid,handles.vidframe))
axes(handles.axes2)
plot(x,0)
xlim([xmin xmax])
ylim([0 1])
vline(handles.vidframe+(str2num(handles.selected(9:end-4))-1)*1000,'k')
vline(handles.grooming(:,1),'g')
vline(handles.grooming(:,2),'r')


% --- Executes on button press in Savebutton.
function Savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to Savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
grooming = handles.grooming;
save('Groomingframes.mat','grooming')

% --- Executes on button press in DeleteStart.
function DeleteStart_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if length(handles.grooming(:,1)) == 1
    handles.grooming(1,:) = 0;
else
    handles.grooming(end,:) = [];
end
guidata(hObject,handles)
updateAxes(handles)

% --- Executes on button press in DeleteStop.
function DeleteStop_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.grooming(end,2) = 0;
guidata(hObject,handles)
updateAxes(handles)
