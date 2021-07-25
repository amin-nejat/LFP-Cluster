function varargout = UI_Evaluation(varargin)



gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UI_Evaluation_OpeningFcn, ...
                   'gui_OutputFcn',  @UI_Evaluation_OutputFcn, ...
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


function UI_Evaluation_OpeningFcn(hObject, eventdata, handles, varargin)
controller = Controller.get_instance();
e = Evaluation.get_instance();

rand_index = zeros(length(controller.partitions), length(controller.partitions));
f_measure = zeros(length(controller.partitions), length(controller.partitions));
precision = zeros(length(controller.partitions), length(controller.partitions));
recall = zeros(length(controller.partitions), length(controller.partitions));

for i = 1 : length(controller.partitions)
    rand_index(i, i) = 1;
    f_measure(i, i) = 1;
    precision(i, i) = 1;
    recall(i, i) = 1;
    for j = i + 1 : length(controller.partitions)
        [r, f, p, re] = e.clusterers_similarity(controller.partitions{i}, controller.partitions{j});
        
        rand_index(i, j) = r;
        rand_index(j, i) = r;
        
        f_measure(i, j) = f;
        f_measure(j, i) = f;
        
        precision(i, j) = p;
        precision(j, i) = p;
        
        recall(i, j) = re;
        recall(j, i) = re;
    end
end


handles.output = hObject;
guidata(hObject, handles);

handles.Rand_Index.set('Data', rand_index);
handles.F_Measure.set('Data', f_measure);
handles.Precision.set('Data', precision);
handles.Recall.set('Data', recall);

function varargout = UI_Evaluation_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;
