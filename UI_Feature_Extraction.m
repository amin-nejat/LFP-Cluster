function varargout = UI_Feature_Extraction(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @UI_Feature_Extraction_OpeningFcn, ...
    'gui_OutputFcn',  @UI_Feature_Extraction_OutputFcn, ...
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

function UI_Feature_Extraction_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = UI_Feature_Extraction_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function Sink_Source_Callback(hObject, eventdata, handles)

function Orientation_Nonselective_Callback(hObject, eventdata, handles)

function Power_Beta_Callback(hObject, eventdata, handles)

function Power_LGamma_Callback(hObject, eventdata, handles)

function Power_HGamma_Callback(hObject, eventdata, handles)

function Orientation_Selectivity_Callback(hObject, eventdata, handles)
Sample_Number_Callback(handles.Sample_Number, eventdata, handles);
return

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Sample_Number_Callback(hObject, eventdata, handles)
if handles.Orientation_Selectivity.Value == 1 || handles.Sample_Number.Value == 1
    return
end

c = Controller.get_instance();

if handles.Orientation_Selectivity.Value-1 == 3
    [signal, local_optima, beta, lgamma, hgamma, main_cross, sink_source, wavelet, peaks] = c.all_features('Orien', handles.Sample_Number.Value-1);
    
    if sink_source == 1
        handles.Sink_Source_Text.String = 'Source';
    else
        handles.Sink_Source_Text.String = 'Sink';
    end
    
    axes(handles.Wavelet_Plot);
    plot(1: length(wavelet), wavelet, 'color', 'blue');
    hold on;
    
    axes(handles.Main_Cross_Plot);
    normalized = main_cross{3}; midvolts_indices = main_cross{2}; main_cross_vec = main_cross{1}; time35 = main_cross_vec(2); time65 = main_cross_vec(3);
    plot(1: length(normalized), normalized,'.-', midvolts_indices, normalized(midvolts_indices),'ro', time35, normalized(time35), 'c*', time65, normalized(time65), 'c*','linewidth',2);
    hold on;
    
    axes(handles.Peaks_Plot);
    plot(1: length(smooth(signal, c.smoothing)),smooth(signal, c.smoothing),'.-',local_optima{1},local_optima{2},'ro','linewidth',2);
    hold on;
    
    axes(handles.Signal_Plot);
    plot(1: length(signal), signal, 'color', 'blue');
    hold on;
    
    axes(handles.Beta_Plot);
    plot(1: length(beta), beta, 'color', 'blue');
    hold on;
    
    axes(handles.LGamma_Plot);
    plot(1: length(lgamma), lgamma, 'color', 'blue');
    hold on;
    
    axes(handles.HGamma_Plot);
    plot(1: length(hgamma), hgamma, 'color', 'blue');
    hold on;
    
    [signal, local_optima, beta, lgamma, hgamma, main_cross, sink_source, wavelet, peaks] = c.all_features('Nonorien', handles.Sample_Number.Value-1);
    
    if sink_source == 1
        handles.Sink_Source_Text.String = 'Source';
    else
        handles.Sink_Source_Text.String = 'Sink';
    end
    
    axes(handles.Wavelet_Plot);
    plot(1: length(wavelet), wavelet, 'color', 'red');
    hold off;
    
    axes(handles.Main_Cross_Plot);
    normalized = main_cross{3}; midvolts_indices = main_cross{2}; main_cross_vec = main_cross{1}; time35 = main_cross_vec(2); time65 = main_cross_vec(3);
    plot(1: length(normalized), normalized,'.-', midvolts_indices, normalized(midvolts_indices),'ro', time35, normalized(time35), 'c*', time65, normalized(time65), 'c*','linewidth',2);
    hold off;
    
    axes(handles.Peaks_Plot);
    plot(1: length(smooth(signal, c.smoothing)),smooth(signal, c.smoothing),'.-',local_optima{1},local_optima{2},'ro','linewidth',2);
    hold off;
    
    axes(handles.Signal_Plot);
    plot(1: length(signal), signal, 'color', 'red');
    hold off;
    
    axes(handles.Beta_Plot);
    plot(1: length(beta), beta, 'color', 'red');
    hold off;
    
    axes(handles.LGamma_Plot);
    plot(1: length(lgamma), lgamma, 'color', 'red');
    hold off;
    
    axes(handles.HGamma_Plot);
    plot(1: length(hgamma), hgamma, 'color', 'red');
    hold off;
    
else
    states = {'Orien', 'Nonorien', 'Both', 'Mean'};
    [signal, local_optima, beta, lgamma, hgamma, main_cross, sink_source, wavelet, peaks] = c.all_features(states{handles.Orientation_Selectivity.Value-1}, handles.Sample_Number.Value-1);
    
    if sink_source == 1
        handles.Sink_Source_Text.String = 'Source';
    else
        handles.Sink_Source_Text.String = 'Sink';
    end
    
    axes(handles.Wavelet_Plot);
    plot(1: length(wavelet), wavelet, 'color', 'blue');

    axes(handles.Main_Cross_Plot);
    normalized = main_cross{3}; midvolts_indices = main_cross{2}; main_cross_vec = main_cross{1}; time35 = main_cross_vec(2); time65 = main_cross_vec(3);
    plot(1: length(normalized), normalized,'.-', midvolts_indices, normalized(midvolts_indices),'ro', time35, normalized(time35), 'c*', time65, normalized(time65), 'c*','linewidth',2);
    
    axes(handles.Peaks_Plot);
    plot(1: length(smooth(signal, c.smoothing)),smooth(signal, c.smoothing),'.-',local_optima{1},local_optima{2},'ro','linewidth',2);
    
    axes(handles.Signal_Plot);
    plot(1: length(signal), signal, 'color', 'blue');
    
    axes(handles.Beta_Plot);
    plot(1: length(beta), beta, 'color', 'blue');
    
    axes(handles.LGamma_Plot);
    plot(1: length(lgamma), lgamma, 'color', 'blue');
    
    axes(handles.HGamma_Plot);
    plot(1: length(hgamma), hgamma, 'color', 'blue');
end

function Sample_Number_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function figure1_CloseRequestFcn(hObject, eventdata, handles)
if handles.Orientation_Selectivity.Value == 1
    return
end

controller = Controller.get_instance();
feature_names = {};

states = {'Orien', 'Nonorien', 'Both', 'Mean'};


if handles.Peaks.Value
    feature_names{length(feature_names)+1} = 'Peaks';
end
if handles.Integral.Value
    feature_names{length(feature_names)+1} = 'Integral';
end
if handles.Main_Cross.Value
    feature_names{length(feature_names)+1} = 'Main_Cross';
end
if handles.Sink_Source.Value
    feature_names{length(feature_names)+1} = 'Sink_Source';
end
if handles.Power_HGamma.Value
    feature_names{length(feature_names)+1} = 'Power_HGamma';
end
if handles.Power_LGamma.Value
    feature_names{length(feature_names)+1} = 'Power_LGamma';
end
if handles.Power_Beta.Value
    feature_names{length(feature_names)+1} = 'Power_Beta';
end
if handles.Wavelet.Value
    feature_names{length(feature_names)+1} = 'Wavelet';
end
controller.construct_fvectors(states{handles.Orientation_Selectivity.Value-1}, feature_names);
delete(hObject);


function Orientation_Selectivity_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Peaks_Callback(hObject, eventdata, handles)




function Smooth_Value_Callback(hObject, eventdata, handles)



function Smooth_Value_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Integral_Callback(hObject, eventdata, handles)



function Main_Cross_Callback(hObject, eventdata, handles)


function Sink_Source_Text_CreateFcn(hObject, eventdata, handles)
[x, map] = imread('Images/Sink_Source.png');
I2 = imresize(x, [200, 600]);
set(hObject, 'cdata', I2);


function Integral_Plot_CreateFcn(hObject, eventdata, handles)
[x, map] = imread('Images/Integral.png');
I2 = imresize(x, [200, 600]);
axes(hObject);
imshow(I2);
set(hObject, 'Tag', 'Integral_Plot');

function Sink_Source_Text_Callback(hObject, eventdata, handles)

function Wavelet_Callback(hObject, eventdata, handles)
