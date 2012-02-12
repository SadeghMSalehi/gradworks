function varargout = image3d(varargin)
% IMAGE3D MATLAB code for image3d.fig
%      IMAGE3D, by itself, creates a new IMAGE3D or raises the existing
%      singleton*.
%
%      H = IMAGE3D returns the handle to a new IMAGE3D or the handle to
%      the existing singleton*.
%
%      IMAGE3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGE3D.M with the given input arguments.
%
%      IMAGE3D('Property','Value',...) creates a new IMAGE3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before image3d_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to image3d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help image3d

% Last Modified by GUIDE v2.5 30-Jan-2012 12:12:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @image3d_OpeningFcn, ...
                   'gui_OutputFcn',  @image3d_OutputFcn, ...
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



% --- Executes just before image3d is made visible.
function image3d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to image3d (see VARARGIN)

% Choose default command line output for image3d
handles.output = hObject;
handles.imgProc = @(x)(x);

if (length(varargin) == 1)
    handles.image = varargin{1};
    handles.szImage = size(handles.image);
    set(handles.slider1,'Min',1);
    set(handles.slider1,'Value',handles.szImage(3)/2);
    guidata(hObject, handles);
elseif (length(varargin) == 4) 
    % Xm, Xv, Xw, and szImage
    handles.ShapeInfo = varargin;
    handles.szImage = varargin{4};
    Xw = varargin{3};
    C0 = zeros(size(Xw));
    C0(1) = 1;
    set(handles.edit1, 'String', mat2str(C0(:)));
end

set(handles.slider1,'Value',handles.szImage(3)/2);
update_Figure(hObject,handles);
% Update handles structure
guidata(hObject, handles);


% UIWAIT makes image3d wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = image3d_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_Figure(hObject, handles);

guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in radiobutton1.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
set(handles.slider1,'Value',handles.szImage(3)/2);
update_Figure(hObject,handles);
guidata(hObject, handles);

% --- Executes on button press in radiobutton2.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
set(handles.slider1,'Value',handles.szImage(1)/2);
update_Figure(hObject,handles);
guidata(hObject, handles);

% --- Executes on button press in radiobutton3.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
set(handles.slider1,'Value',handles.szImage(2)/2);
update_Figure(hObject,handles);
guidata(hObject, handles);


   
function update_Figure(hObject, handles)
    coeff = get(handles.edit1, 'String');
    C = eval(coeff);
    if (existsfield(handles,'ShapeInfo'))
        handles.image = reconPCA3d(C,handles.ShapeInfo{1},handles.ShapeInfo{2},handles.ShapeInfo{3},handles.szImage);
    end
    
    imgScale = [min(handles.image(:)) max(handles.image(:))];
    pos = uint8(get(handles.slider1, 'Value'));
    axis = 3;
    if (get(handles.radiobutton4, 'Value') > 0)
        set(handles.slider1,'Max',uint8(handles.szImage(3)));
        pos = min(max(1,pos), handles.szImage(3));
        outImage = (reshape(handles.image(:,:,pos), [handles.szImage(1) handles.szImage(2)]));
        axis = 3;
    elseif get(handles.radiobutton5, 'Value') > 0
            set(handles.slider1,'Max',uint8(handles.szImage(1)));
            pos = min(max(1,pos), handles.szImage(1));
            outImage = (reshape(handles.image(pos,:,:), [handles.szImage(2) handles.szImage(3)]));
            axis = 1;
        elseif get(handles.radiobutton6, 'Value') > 0
            set(handles.slider1,'Max',uint8(handles.szImage(2)));
            pos = min(max(1,pos), handles.szImage(2));
            outImage = (reshape(handles.image(:,pos,:), [handles.szImage(1) handles.szImage(3)]));
            axis = 2;
    end
    set(handles.slider1, 'SliderStep', [1 1/handles.szImage(axis)]);
    set(handles.text1, 'String', sprintf('%03d', pos));
    if (get(handles.checkbox1, 'Value') > 0)
        outImage(outImage>0)  = 1;
        outImage(outImage<=0) = 0;
        outImage = 1-outImage;
        imgScale = [0 1];
    end
    imagesc(outImage, imgScale);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
    update_Figure(hObject,handles);
    guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
guidata(hObject, handles);
