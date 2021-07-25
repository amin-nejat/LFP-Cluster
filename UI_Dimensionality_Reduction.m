function varargout = UI_Dimensionality_Reduction(varargin)



gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UI_Dimensionality_Reduction_OpeningFcn, ...
                   'gui_OutputFcn',  @UI_Dimensionality_Reduction_OutputFcn, ...
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


function UI_Dimensionality_Reduction_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

guidata(hObject, handles);



function varargout = UI_Dimensionality_Reduction_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


function Reduce_Dimension_Callback(hObject, eventdata, handles)
dr = Dimensionality_Reduction.get_instance();
c = Clustering.get_instance();
dr.new_dim = str2num(get(handles.New_Dimension,'String'));
c.X = dr.pca_reduction();


function New_Dimension_Callback(hObject, eventdata, handles)



function New_Dimension_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Visualize_pca_Callback(hObject, eventdata, handles)
dr = Dimensionality_Reduction.get_instance();
dr.pca_visualize();