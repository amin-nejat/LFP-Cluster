clear all;
close all;
clc;

data_folderpath = 'C:\Users\Amin\Desktop\University\Term9\Neuroscience\Project\Data\';
controller = Controller.get_instance();
controller.set_folderpath(data_folderpath);

% controller.fit_waveforms()

UI();