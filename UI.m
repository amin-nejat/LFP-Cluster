function varargout = UI(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @UI_OpeningFcn, ...
    'gui_OutputFcn',  @UI_OutputFcn, ...
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


function UI_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);


function varargout = UI_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


function Feature_Extraction_Callback(hObject, eventdata, handles)
UI_Feature_Extraction();

function Dimensionality_Reduction_Callback(hObject, eventdata, handles)
UI_Dimensionality_Reduction();

function Clustering_Callback(hObject, eventdata, handles)
UI_Clustering();


function Evaluation_Callback(hObject, eventdata, handles)
UI_Evaluation();


function Feature_Extraction_CreateFcn(hObject, eventdata, handles)
[x,map]=imread('Images/Feature_Extraction.png');
pos = get(hObject, 'position');
I2=imresize(x, [200, 200]);
set(hObject, 'cdata', I2);


function Dimensionality_Reduction_CreateFcn(hObject, eventdata, handles)
[x,map]=imread('Images/Dimensionality_Reduction.png');
pos = get(hObject, 'position');
I2=imresize(x, [200, 200]);
set(hObject, 'cdata', I2);


function Clustering_CreateFcn(hObject, eventdata, handles)
[x,map]=imread('Images/Clustering.png');
pos = get(hObject, 'position');
I2=imresize(x, [200, 200]);
set(hObject, 'cdata', I2);


function Evaluation_CreateFcn(hObject, eventdata, handles)
[x,map]=imread('Images/Evaluation.png');
pos = get(hObject, 'position');
I2=imresize(x, [200, 200]);
set(hObject, 'cdata', I2);
