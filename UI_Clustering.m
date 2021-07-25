function varargout = UI_Clustering(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UI_Clustering_OpeningFcn, ...
                   'gui_OutputFcn',  @UI_Clustering_OutputFcn, ...
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


function UI_Clustering_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
guidata(hObject, handles);



function varargout = UI_Clustering_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


function Choose_K_Kmeans_Callback(hObject, eventdata, handles)
    c = Dimensionality_Reduction.get_instance();
    c.kmeans_getk();

function Choose_K_svd_Callback(hObject, eventdata, handles)
    c = Dimensionality_Reduction.get_instance();
    c.svd_getk();

function Cluster_Callback(hObject, eventdata, handles)
c = Clustering.get_instance();
c.k = str2num(get(handles.K, 'String'));
method = get(get(handles.Clustering_Method,'SelectedObject'), 'Tag');
if strcmp(method, 'method_kmeans')
    c.labels = c.kmeans_clustering();
    c.visualize_clustering();
elseif strcmp(method, 'method_gmm')
    c.labels = c.gmm_clustering();
    c.visualize_clustering();
elseif strcmp(method, 'method_bce')
    idx_gmm = c.gmm_clustering();
    idx_kmeans = c.kmeans_clustering();
    c.labels = c.bce_clustering([idx_gmm, idx_kmeans]);
    c.visualize_clustering();
end
controller = Controller.get_instance();
controller.partitions{length(controller.partitions) + 1} = c.labels;

function K_Callback(hObject, eventdata, handles)


function K_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
