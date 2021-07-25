classdef Controller < handle
    properties
        files
        ori_folderpath
        smoothing = 500
        signal_length = 1000
        partitions = {[1;1;1;2;2;2;3;3;3;3;3;3;4;4;5;5;5;5;6;6;6;7;7;7;7;7;7;7;7;8;8;8;9;9;9;9;9;9;10;10;10;10]}
        orientation_selectivity = 1
    end
    
    methods (Access = private)
        function obj = Controller
            data_folderpath = 'C:\Users\Amin\Desktop\University\Term9\Neuroscience\Project\Data\';
            obj.ori_folderpath = strcat(data_folderpath, 'Orientation_data\Orientation\');
            obj.files = dir([obj.ori_folderpath '*.mat']);
        end
    end
    
    methods(Static)
        function r = get_instance
            persistent instance;
            if isempty(instance) || ~isvalid(instance)
                instance = Controller();
            end
            r = instance;
        end
    end
    
    methods
        function set_folderpath(obj, data_folderpath)
            obj.ori_folderpath = strcat(data_folderpath, 'Orientation_data\Orientation\');
            obj.files = dir([obj.ori_folderpath '*.mat']);
        end
        
        function[signal, local_optima, beta, lgamma, hgamma, main_cross, sink_source, wavelet, peaks] = all_features(obj, state, index)
            
            filepath = strcat(obj.ori_folderpath, obj.files(index).name);
            load(filepath);
            
            feature_extraction = Feature_Extraction.get_instance();
            [orientations, mean_orien, mean_nonorien, mean_total] = feature_extraction.get_orientations(neuronlfp);
            
            if strcmp(state, 'Orien')
                signal = mean_orien;
            elseif strcmp (state, 'Nonorien')
                signal = mean_nonorien;
            elseif strcmp (state, 'Mean')
                signal = mean_total;
            end
            
            [max_indices, max_loc] = feature_extraction.peakfinder(smooth(signal, obj.smoothing), 1);
            [min_indices, min_loc] = feature_extraction.peakfinder(smooth(signal, obj.smoothing), -1);
            
            peaks = [length(min_indices), length(max_indices)];
            
            sink_source = feature_extraction.sink_source(min_indices, min_loc, max_indices, max_loc, smooth(signal, obj.smoothing));
            
            local_optima_indices = [max_indices; min_indices];
            local_optima_loc = [max_loc; min_loc];
            local_optima = {local_optima_indices, local_optima_loc};
            
            beta = real(feature_extraction.lfp_fft_beta(signal));
            lgamma = real(feature_extraction.lfp_fft_lgamma(signal));
            hgamma = real(feature_extraction.lfp_fft_hgamma(signal));
            
            [shape, time35, time65, midvolts, midvolts_indices, normalized] = feature_extraction.main_cross(smooth(signal, obj.smoothing));
            main_cross = {[shape, time35, time65, midvolts], midvolts_indices, normalized};
            
            wavelet = feature_extraction.wavelet_pulse(signal);
        end
        
        function fit_waveforms(obj)
            basis_mat = [];
            for i = 1 : length(obj.files)
                filepath = strcat(obj.ori_folderpath, obj.files(i).name);
                load(filepath);
                wave = transpose(mean(neuronwaf{1,1}, 2));
                
                if isempty(basis_mat)
                    basis_mat = wave;
                elseif any(isnan(wave)) == 0
                    basis_mat = [basis_mat ; wave];
                end
                
                
            end
            size(mean(neuronlfp{1}))
            feature_extraction = Feature_Extraction.get_instance();
            feature_extraction.fit_waveforms(neuronlfp{1}, basis_mat);
            
        end
        
        function [signal, fvector] = construct_fvector(obj, state, feature_names, index)
            feature_extraction = Feature_Extraction.get_instance();
            [signal, local_optima, beta, lgamma, hgamma, main_cross, sink_source, wavelet, peaks] = obj.all_features(state, index);
            tmp = [];
            if any(strcmp(feature_names, 'Sink_Source'))
                tmp = [tmp, sink_source];
            end
            if any(strcmp(feature_names, 'Peaks'))
                tmp = [tmp, peaks];
            end
            if any(strcmp(feature_names, 'Integral'))
                tmp = [tmp, feature_extraction.integral_pulse(signal)];
            end
            if any(strcmp(feature_names, 'Main_Cross'))
                tmp = [tmp, main_cross{1}];
            end
            if any(strcmp(feature_names, 'Power_Beta'))
                tmp = [tmp, real(beta)];
            end
            if any(strcmp(feature_names, 'Power_LGamma'))
                tmp = [tmp, real(lgamma)];
            end
            if any(strcmp(feature_names, 'Power_HGamma'))
                tmp = [tmp, real(hgamma)];
            end
            if any(strcmp(feature_names, 'Wavelet'))
                tmp = [tmp, wavelet];
            end
            fvector = tmp;
        end
        
        function construct_fvectors(obj, state, feature_names)
            feature_extraction = Feature_Extraction.get_instance();
            feature_extraction.X = [];
            
            c = Clustering.get_instance();
            c.signals_orien = []; c.signals_nonorien = []; c.signals = [];
            
            dr = Dimensionality_Reduction.get_instance();
            
            for i = 1 : length(obj.files)
                
                if strcmp(state, 'Both')
                    [signal1, fvector1] = obj.construct_fvector('Orien', feature_names, i);
                    [signal2, fvector2] = obj.construct_fvector('Nonorien', feature_names, i);
                    
                    c.signals_orien = [c.signals_orien; signal1];
                    c.signals_nonorien = [c.signals_nonorien; signal2];
                    
                    tmp = [fvector1, fvector2];
                else
                    [signal, tmp] = obj.construct_fvector(state, feature_names, i);
                    c.signals = [c.signals; signal];
                end
                
                

                feature_extraction.X = [feature_extraction.X; tmp];
            end
            
            dr.X = feature_extraction.X;
            c.X = feature_extraction.X;
        end
    end
end