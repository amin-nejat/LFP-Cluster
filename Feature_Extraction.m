classdef Feature_Extraction < handle
    properties
        X
    end
    
    methods (Access = private)
        function obj = Feature_Extraction
        end
    end
    
    methods(Static)
        function r = get_instance
            persistent instance;
            if isempty(instance)
                instance = Feature_Extraction;
            end
            r = instance;
        end
    end
    
    methods
        function peak_feature(obj, x0)
            [max_indices, max_loc] = peakfinder(smooth(x0, 10), 1);
            [min_indices, min_loc] = peakfinder(smooth(x0, 10), -1);
        end
        
        function result = sink_source(obj, min_indices, min_loc, max_indices, max_loc, signal)
           variance = sqrt(var(signal));
           m = mean(signal);
           
           min_loc = abs(m - min_loc);
           max_loc = abs(m - max_loc);
           
           imin = find(min_loc > variance, 1);
           imax = find(max_loc > variance, 1);
           
           result = 2*(min_indices(imin) > max_indices(imax))-1;
           if isempty(result)
               if isempty(max_indices)
                   result = -1;
               elseif isempty(min_indices)
                   result = 1;
               else
                   result = 2*(min_indices(1) > max_indices(1)) - 1;
               end
           end
        end
        
        function varargout = peakfinder(obj, x0, extrema)
            sel = (max(x0)-min(x0))/4;
            includeEndpoints = false;
            interpolate = false;
            thresh = [];
            
            x0 = extrema*x0(:); % Make it so we are finding maxima regardless
            thresh = thresh*extrema; % Adjust threshold according to extrema.
            dx0 = diff(x0); % Find derivative
            dx0(dx0 == 0) = -eps; % This is so we find the first of repeated values
            ind = find(dx0(1:end-1).*dx0(2:end) < 0)+1; % Find where the derivative changes sign
            len0 = numel(x0); s = size(x0); flipData =  s(1) < s(2);
            
            
            % Include endpoints in potential peaks and valleys as desired
            if includeEndpoints
                x = [x0(1);x0(ind);x0(end)];
                ind = [1;ind;len0];
                minMag = min(x);
                leftMin = minMag;
            else
                x = x0(ind);
                minMag = min(x);
                leftMin = min(x(1), x0(1));
            end
            
            % x only has the peaks, valleys, and possibly endpoints
            len = numel(x);
            
            if len > 2 % Function with peaks and valleys
                % Set initial parameters for loop
                tempMag = minMag;
                foundPeak = false;
                
                if includeEndpoints
                    % Deal with first point a little differently since tacked it on
                    % Calculate the sign of the derivative since we tacked the first
                    %  point on it does not neccessarily alternate like the rest.
                    signDx = sign(diff(x(1:3)));
                    if signDx(1) <= 0 % The first point is larger or equal to the second
                        if signDx(1) == signDx(2) % Want alternating signs
                            x(2) = [];
                            ind(2) = [];
                            len = len-1;
                        end
                    else % First point is smaller than the second
                        if signDx(1) == signDx(2) % Want alternating signs
                            x(1) = [];
                            ind(1) = [];
                            len = len-1;
                        end
                    end
                end
                
                % Skip the first point if it is smaller so we always start on a
                %   maxima
                if x(1) >= x(2)
                    ii = 0;
                else
                    ii = 1;
                end
                
                % Preallocate max number of maxima
                maxPeaks = ceil(len/2);
                peakLoc = zeros(maxPeaks,1);
                peakMag = zeros(maxPeaks,1);
                cInd = 1;
                % Loop through extrema which should be peaks and then valleys
                while ii < len
                    ii = ii+1; % This is a peak
                    % Reset peak finding if we had a peak and the next peak is bigger
                    %   than the last or the left min was small enough to reset.
                    if foundPeak
                        tempMag = minMag;
                        foundPeak = false;
                    end
                    
                    % Found new peak that was lager than temp mag and selectivity larger
                    %   than the minimum to its left.
                    if x(ii) > tempMag && x(ii) > leftMin + sel
                        tempLoc = ii;
                        tempMag = x(ii);
                    end
                    
                    % Make sure we don't iterate past the length of our vector
                    if ii == len
                        break; % We assign the last point differently out of the loop
                    end
                    
                    ii = ii+1; % Move onto the valley
                    % Come down at least sel from peak
                    if ~foundPeak && tempMag > sel + x(ii)
                        foundPeak = true; % We have found a peak
                        leftMin = x(ii);
                        peakLoc(cInd) = tempLoc; % Add peak to index
                        peakMag(cInd) = tempMag;
                        cInd = cInd+1;
                    elseif x(ii) < leftMin % New left minima
                        leftMin = x(ii);
                    end
                end
                
                % Check end point
                if includeEndpoints
                    if x(end) > tempMag && x(end) > leftMin + sel
                        peakLoc(cInd) = len;
                        peakMag(cInd) = x(end);
                        cInd = cInd + 1;
                    elseif ~foundPeak && tempMag > minMag % Check if we still need to add the last point
                        peakLoc(cInd) = tempLoc;
                        peakMag(cInd) = tempMag;
                        cInd = cInd + 1;
                    end
                elseif ~foundPeak
                    if x(end) > tempMag && x(end) > leftMin + sel
                        peakLoc(cInd) = len;
                        peakMag(cInd) = x(end);
                        cInd = cInd + 1;
                    elseif tempMag > min(x0(end), x(end)) + sel
                        peakLoc(cInd) = tempLoc;
                        peakMag(cInd) = tempMag;
                        cInd = cInd + 1;
                    end
                end
                
                % Create output
                if cInd > 1
                    peakInds = ind(peakLoc(1:cInd-1));
                    peakMags = peakMag(1:cInd-1);
                else
                    peakInds = [];
                    peakMags = [];
                end
            else % This is a monotone function where an endpoint is the only peak
                [peakMags,xInd] = max(x);
                if includeEndpoints && peakMags > minMag + sel
                    peakInds = ind(xInd);
                else
                    peakMags = [];
                    peakInds = [];
                end
            end
            
            % Apply threshold value.  Since always finding maxima it will always be
            %   larger than the thresh.
            if ~isempty(thresh)
                m = peakMags>thresh;
                peakInds = peakInds(m);
                peakMags = peakMags(m);
            end
            
            if interpolate && ~isempty(peakMags)
                middleMask = (peakInds > 1) & (peakInds < len0);
                noEnds = peakInds(middleMask);
                
                magDiff = x0(noEnds + 1) - x0(noEnds - 1);
                magSum = x0(noEnds - 1) + x0(noEnds + 1)  - 2 * x0(noEnds);
                magRatio = magDiff ./ magSum;
                
                peakInds(middleMask) = peakInds(middleMask) - magRatio/2;
                peakMags(middleMask) = peakMags(middleMask) - magRatio .* magDiff/8;
            end
            
            % Rotate data if needed
            if flipData
                peakMags = peakMags.';
                peakInds = peakInds.';
            end
            
            % Change sign of data if was finding minima
            if extrema < 0
                peakMags = -peakMags;
                x0 = -x0;
            end
            
            % Plot if no output desired
            if nargout == 0
                if isempty(peakInds)
                    disp('No significant peaks found')
                else
                    figure;
                    plot(1:len0,x0,'.-',peakInds,peakMags,'ro','linewidth',2);
                end
            else
                varargout = {peakInds,peakMags};
            end
        end
        
        
        function Beta = lfp_fft_beta(obj, signal)
            
            fa=12; fb=30;
            FT=fft(signal);
            Beta=zeros(1,floor(0.2*fb)-floor(0.2*fa)+1);
            
            for h=floor(0.2*fa):floor(0.2*fb)
                Beta(h-floor(0.2*fa)+1)=FT(h);
            end
        end
        
        
        
        
        function LGamma = lfp_fft_lgamma(obj, signal)
            
            fa=30; fb=90;
            FT=fft(signal);
            
            LGamma=zeros(1,floor(0.2*fb)-floor(0.2*fa)+1);
            
            for h=floor(0.2*fa):floor(0.2*fb)
                LGamma(h-floor(0.2*fa)+1)=FT(h);
            end
            
        end
        
        
        function HGamma = lfp_fft_hgamma(obj, signal)
            
            fa=90; fb=200;
            FT=fft(signal);
            HGamma=zeros(1,floor(0.2*fb)-floor(0.2*fa)+1);
            
            for h=floor(0.2*fa):floor(0.2*fb)
                HGamma(h-floor(0.2*fa)+1)=FT(h);
            end
            
        end
        
        function wavelet = wavelet_pulse(obj, signal)
           wavelet = dwt(signal, 'db5'); 
        end
        
        function Integral = integral_pulse(obj,signal)  %According to the paper with DOI:  10.1152/jn.00278.2014
            Integral=sum(signal);
        end
        
        
        
        function Range= range_pulse(obj,signal)  %According to the paper with DOI:  10.1152/jn.00278.2014
            Range=max(signal)-min(signal);
        end
        
        
        %        function TimeFirstPeak=time_first_peak(obj,signal)%According to the paper with DOI:  10.1152/jn.00278.2014
        %
        %            %%in code ro khodet ba peakfinder fek konam bezani behtare
        %
        %        end
        %
        
        
        
        function [shape1, time35, time65, MidVolts, midvolts_indices, signal]= main_cross(obj, signal)
            
%             signal = smooth(signal,60);  %smoothing
            signal = signal/max(abs(signal(50:600)));% normalize, ignoring the edges
            
            xMax = find(signal == max(signal(50:600)), 1, 'first');
            xMin = find(signal == min(signal(50:600)), 1, 'first');
            
            if xMax>xMin
                shape1 = -1;
            else
                shape1 = 1;
            end
            
            if shape1 == -1
                
                for i = xMin : xMax
                    if(signal(i) - signal(xMin) > 0.35*(signal(xMax) - signal(xMin)))
                        time35 = i;
                        break;
                    end
                end
                
                for i = xMax : -1 : xMin
                    if(signal(i) - signal(xMin) < 0.65*(signal(xMax) - signal(xMin)))
                        time65 = i;
                        break;
                    end
                end
            end
            
            if shape1 == 1
                
                for i = xMax : xMin
                    if (signal(i) - signal(xMin) < 0.65*(signal(xMax) - signal(xMin)))
                        time35 = i;
                        break;
                    end
                end
                
                for i = xMin : -1 : xMax
                    if(signal(i) - signal(xMin) > 0.35*(signal(xMax) - signal(xMin)))
                        time65 = i;
                        break;
                    end
                end
                
                
            end
            
            times = linspace(time35,time65,10);
            MidVolts = zeros(1,10);
            
            for l = 1 : 10
                MidVolts(l) = (signal(ceil(times(l))) - signal(floor(times(l)))) * (times(l) - floor(times(l))) + signal(floor(times(l)));   %f=(fb-fb)(x-xa)+fa
            end
            
            midvolts_indices = ceil(times);
            
        end
        
        
        
        function[distance, x, y] = signals_distance(obj, f, g)
            
            if length(f) ~= length(g)
                disp('Lengths must be the same');
                return
            end
            
            a = sum(f);
            b = sum(g);
            c = dot(f, f);
            d = dot(f, g);
            k = length(f);
            
            x = (a*d-b*c)/(a*b-d*k);
            y = (a*b-d*k)/(a*a-c*k);
            
            distance = dot((x+f)*y-g, (x+f)*y-g);
            
            
            %     figure(1);
            %     plot(linspace(0, 200, length(f)), (x+f)*y, 'color', 'red');
            %     hold on;
            %     plot(linspace(0, 200, length(f)), g, 'color', 'blue');
            %     hold on;
            %     plot(linspace(0, 200, length(f)), f, 'color', 'yellow');
            %     hold on;
            %
            %     syms x y z t;
            %     r = solve([a*y*y+k*y*(x*y+z*t)+y*t*b == 0, b*t*t+k*t*(x*y+z*t)+y*t*a == 0, y*c+2*x*y*a+k*x*(x*y+z*t)+t*e+x*t*b+z*t*a == 0, t*d+2*z*t*b+k*z*(x*y+z*t)+y*e+z*y*a+x*y*b == 0], [x,y,z,t])
            
        end
        
        
        function[orientations, mean_orien, mean_nonorien, mean_total] = get_orientations(obj, neuronlfp)
            
            sim = zeros(16, 16);
            indices = zeros(1, 16);
            coeff_x = zeros(16, 16);
            coeff_y = eye(16, 16);
            for m = 1 : 16
                for n = m+1 : 16
                    [sim(m, n), coeff_y(m, n), coeff_y(m, n)] = obj.signals_distance(mean(neuronlfp{m}(:,1001:2000), 1), mean(neuronlfp{n}(:,1001:2000), 1));
                    sim(n, m) = sim(m, n);
                end
                [m, index] = max(sim(m,:));
                indices(index) = indices(index)+1;
            end
            [m, index] = sort(sum(sim), 'descend');
            orientations = zeros(1,16);
            for i = 1 : 16
                if indices(index(i)) > 0
                    orientations(i) = index(i);
                else
                    break
                end
            end
            orientations = orientations(orientations>0);
            
            mean_orien = mean(neuronlfp{1}(:,1001:2000), 1);
            mean_nonorien = mean(neuronlfp{2}(:,1001:2000), 1);
            for i = 1 : 16
                if any(i == orientations)
                    mean_orien = mean_orien + coeff_y(1, i)*(coeff_x(1, i)+mean(neuronlfp{i}(:,1001:2000), 1));
                    %plot(linspace(0, 200, length(mean(neuronlfp{i}(:,1001:2000), 1))), coeff_y(1, i)*(coeff_x(1, i)+mean(neuronlfp{i}(:,1001:2000), 1)), 'color', 'red');
                else
                    mean_nonorien = mean_nonorien + coeff_y(1, i)*(coeff_x(1, i)+mean(neuronlfp{i}(:,1001:2000), 1));
                    %plot(linspace(0, 200, length(mean(neuronlfp{i}(:,1001:2000), 1))), coeff_y(1, i)*(coeff_x(1, i)+mean(neuronlfp{i}(:,1001:2000), 1)), 'color', 'blue');
                end
                %         hold on;
            end
            mean_orien = mean_orien/length(orientations);
            mean_nonorien = mean_nonorien/(16-length(orientations));
            mean_total = (mean_nonorien*length(mean_nonorien) + mean_orien*length(mean_orien))/(length(mean_orien) + length(mean_nonorien));
            
            %mean_orien = smooth(mean_orien, 20);
            %mean_nonorien = smooth(mean_nonorien, 20);
            %mean_total = smooth(mean_total, 20);
            
            %             figure;
            %             plot(linspace(0, 200, length(mean_orien)), mean_orien, 'color', 'blue');
            %             hold on;
            %             plot(linspace(0, 200, length(mean_nonorien)), mean_nonorien, 'color', 'red');
            %             hold on;
            %             plot(linspace(0, 200, length(mean_nonorien)), (mean_nonorien*length(mean_nonorien) + mean_orien*length(mean_orien))/(length(mean_orien) + length(mean_nonorien)), 'color', 'yellow');
            %             hold on;
            
        end
        function fit_waveforms(obj, signal, waveforms)
            % x = sym('x', [1, size(basis_mat, 1)]);
            % assume(x, 'real');
            %
            % y = sym('y', [1, size(basis_mat, 1)]);
            % assume(y, 'real');
            %
            % fft_x = real(feval(symengine, 'numeric::fft', x*basis_mat+y));
            basis_mat = waveforms;
            signal = signal(1:size(basis_mat, 2));
            [x_hat, fval_hat] = fminunc(@(x)cost_func(x), 10*rand(1, 2*size(basis_mat, 1)));
            
            x_hat(size(basis_mat, 1):end) = ceil(x_hat(size(basis_mat, 1):end));
            tmp_hat = zeros(size(basis_mat, 1), size(basis_mat, 2));
            for i = 1 : size(basis_mat, 1)
                tmp_hat(i, :) = [zeros(1, x_hat(size(basis_mat, 1) + i)), basis_mat(i, 1:size(basis_mat, 2) - x_hat(size(basis_mat, 1) + i))];
            end
            signal
            figure(1);
%             plot(1:length(signal), signal);
%             hold on;
            result = x_hat(1:size(basis_mat, 1))*tmp_hat;
            plot(1:length(result), result);
            
            function result = cost_func(x)
                x(size(basis_mat, 1):end) = ceil(x(size(basis_mat, 1):end));
                tmp = zeros(size(basis_mat, 1), size(basis_mat, 2));
                for i = 1 : size(basis_mat, 1)
                    tmp(i, :) = [zeros(1, x(size(basis_mat, 1) + i)), basis_mat(i, 1:size(basis_mat, 2) - x(size(basis_mat, 1) + i))];
                end

                fft_x = real(fft(x(1:size(basis_mat, 1))*tmp));
                result = dot(fft_x - fft(signal), fft_x - fft(signal));
            end
        end
        
    end
end