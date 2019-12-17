function varargout = cellReg_2019(varargin)
% CELLREG_2019 MATLAB code for cellReg_2019.fig
%      CELLREG_2019, by itself, creates a new CELLREG_2019 or raises the existing
%      singleton*.
%
%      H = CELLREG_2019 returns the handle to a new CELLREG_2019 or the handle to
%      the existing singleton*.
%
%      CELLREG_2019('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLREG_2019.M with the given input arguments.
%
%      CELLREG_2019('Property','Value',...) creates a new CELLREG_2019 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cellReg_2019_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cellReg_2019_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cellReg_2019

% Last Modified by GUIDE v2.5 05-Nov-2018 16:39:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cellReg_2019_OpeningFcn, ...
                   'gui_OutputFcn',  @cellReg_2019_OutputFcn, ...
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


% --- Executes just before cellReg_2019 is made visible.
function cellReg_2019_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cellReg_2019 (see VARARGIN)

% Choose default command line output for cellReg_2019
handles.output = hObject;
handles.Technique = 'one_photon';
handles.micronsperpixel = 2.35;
handles.Alignment = 'Non-rigid';
handles.maxdist = 12;
handles.maxrotations = 30;
handles.rotationsmooth = 2;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cellReg_2019 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cellReg_2019_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in RunButton.
function RunButton_Callback(hObject, eventdata, handles)
% hObject    handle to RunButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LandmarkShift(handles.folderpath);
alignPairwiseSessionsGUI(handles.folderpath, handles.Technique, handles.micronsperpixel, handles.Alignment,handles.maxdist,handles.maxrotations, handles.rotationsmooth)



% --- Executes on button press in LoadData.
function LoadData_Callback(hObject, eventdata, handles)
% hObject    handle to LoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folderpath = uigetdir(pwd,'Select Folder Containing all sessions miniscope variables');
mkdir(folderpath, 'OriginalFiles')
copyfile(folderpath,[folderpath,'/OriginalFiles']);
handles.folderpath = folderpath;
CurrentFolderView_Callback(hObject,eventdata,handles);
guidata(hObject,handles)


% --- Executes on selection change in CurrentFolderView.
function CurrentFolderView_Callback(hObject, eventdata, handles)
% hObject    handle to CurrentFolderView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.contents = dir([handles.folderpath]);
vlist = {handles.contents(~[handles.contents.isdir]).name};
set(handles.CurrentFolderView,'String',vlist)
guidata(hObject,handles)
% Hints: contents = cellstr(get(hObject,'String')) returns CurrentFolderView contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CurrentFolderView


% --- Executes during object creation, after setting all properties.
function CurrentFolderView_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrentFolderView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OnePhotonButton.
function OnePhotonButton_Callback(hObject, eventdata, handles)
% hObject    handle to OnePhotonButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') == 1
    set(handles.TwoPhotonButton,'Value',0)
    handles.Technique = 'one_photon';
end
guidata(hObject,handles)
    
% Hint: get(hObject,'Value') returns toggle state of OnePhotonButton


% --- Executes on button press in TwoPhotonButton.
function TwoPhotonButton_Callback(hObject, eventdata, handles)
% hObject    handle to TwoPhotonButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') == 1
    set(handles.OnePhotonButton,'Value',0)
    handles.Technique = 'two_photon';    
end
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of TwoPhotonButton



function P_sameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to P_sameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Psame = str2double(get(hObject,'String'));
guidata(hObject,handles)
% Hints: get(hObject,'String') returns contents of P_sameEdit as text
%        str2double(get(hObject,'String')) returns contents of P_sameEdit as a double


% --- Executes during object creation, after setting all properties.
function P_sameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to P_sameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxDistEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MaxDistEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.maxdist = str2double(get(hObject,'String'));
guidata(hObject,handles)
% Hints: get(hObject,'String') returns contents of MaxDistEdit as text
%        str2double(get(hObject,'String')) returns contents of MaxDistEdit as a double


% --- Executes during object creation, after setting all properties.
function MaxDistEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxDistEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TranslationsButton.
function TranslationsButton_Callback(hObject, eventdata, handles)
% hObject    handle to TranslationsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') == 1
    set(handles.TranslationsRotationsButton,'Value',0)
    set(handles.NonRigidButton,'Value',0)
    handles.Technique = 'Translations';   
end
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of TranslationsButton


% --- Executes on button press in TranslationsRotationsButton.
function TranslationsRotationsButton_Callback(hObject, eventdata, handles)
% hObject    handle to TranslationsRotationsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') == 1
    set(handles.TranslationsButton,'Value',0)
    set(handles.NonRigidButton,'Value',0)
    handles.Alignment = 'Translations and Rotations';   
end
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of TranslationsRotationsButton


% --- Executes on button press in NonRigidButton.
function NonRigidButton_Callback(hObject, eventdata, handles)
% hObject    handle to NonRigidButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') == 1
    set(handles.TranslationsRotationsButton,'Value',0)
    set(handles.TranslationsButton,'Value',0)
    handles.Technique = 'Non-rigid';   
end
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of NonRigidButton



function DegRotationsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to DegRotationsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.maxrotations = str2double(get(hObject,'String'));
guidata(hObject,handles)
% Hints: get(hObject,'String') returns contents of DegRotationsEdit as text
%        str2double(get(hObject,'String')) returns contents of DegRotationsEdit as a double


% --- Executes during object creation, after setting all properties.
function DegRotationsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DegRotationsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TransSmoothEdit_Callback(hObject, eventdata, handles)
% hObject    handle to TransSmoothEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.rotationsmooth = str2double(get(hObject,'String'));
guidata(hObject,handles)
% Hints: get(hObject,'String') returns contents of TransSmoothEdit as text
%        str2double(get(hObject,'String')) returns contents of TransSmoothEdit as a double


% --- Executes during object creation, after setting all properties.
function TransSmoothEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TransSmoothEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MicPerPixelEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MicPerPixelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.micronsperpixel = str2double(get(hObject,'String'));
guidata(hObject,handles)
% Hints: get(hObject,'String') returns contents of MicPerPixelEdit as text
%        str2double(get(hObject,'String')) returns contents of MicPerPixelEdit as a double


% --- Executes during object creation, after setting all properties.
function MicPerPixelEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MicPerPixelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
