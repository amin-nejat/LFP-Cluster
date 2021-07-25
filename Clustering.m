classdef Clustering < handle
    properties
        X
        k
        labels
        signals
        signals_orien
        signals_nonorien
        
    end
    
    methods (Access = private)
      function obj = Clustering
      end
    end
    
    methods(Static)
      function r = get_instance
         persistent instance;
         if isempty(instance)
            instance = Clustering(); 
         end
         r = instance;
      end
    end
    
    methods
        
        
        function r = kmeans_clustering(obj)
            [idx, C, sumd] = kmeans(obj.X, obj.k);
            r = idx;
        end
        
        
        function r = gmm_clustering(obj)
            GMM = fitgmdist(obj.X, obj.k, 'RegularizationValue', 0.01, 'Start', 'plus');
            r = GMM.cluster(obj.X);
        end
        
        
        function visualize_clustering(obj)
            k = max(obj.labels);
            
            for i = 1 : k
                figure(i);
                current_pos = 1;
                div_x = ceil(sqrt(sum(obj.labels == i)));
                div_y = ceil(sum(obj.labels == i)/div_x);
                for signal_index = 1 : length(obj.labels)
                    if obj.labels(signal_index) == i
                        subplot(div_x, div_y, current_pos);
                        if isempty(obj.signals)
                            plot(linspace(0, 200, length(obj.signals_orien(signal_index, :))), obj.signals_orien(signal_index, :), 'color', 'red');
                            hold on;
                            
                            plot(linspace(0, 200, length(obj.signals_nonorien(signal_index, :))), obj.signals_nonorien(signal_index, :), 'color', 'blue');
                            hold on;
                        else
                            plot(linspace(0, 200, length(obj.signals(signal_index, :))), obj.signals(signal_index, :), 'color', 'blue');
                            hold on;
                        end
                        current_pos = current_pos + 1;
                    end
                end
            end
        end
        
        
        function r = bce_clustering(obj, base_labels)
            k = max(base_labels(:, 1));
            q = k;
            Palpha = rand(k,1);
            
            for i = 1 : size(base_labels, 2)
                temp = rand(k, q);
                [k,q] = size(temp);
                temp=temp./(sum(temp,2)*ones(1,q));
                Pbeta{i}=temp;
            end
            [phiAll, gammaAll, resultAlpha, resultBeta] = obj.bce_learn(base_labels, Palpha, Pbeta, 0.000001, max(base_labels));
            
            r = zeros(1, size(base_labels, 1));
            for index = 1 : size(base_labels, 1)
                wtheta(:, index) = gammaAll(:, index);
                bb = find(wtheta(:, index) == max(wtheta(:, index)));
                r(index) = bb(1);
            end
        end
        
        
        function [phiAll,gamaAll,resultAlpha,resultBeta] = bce_learn(obj, X, oldAlpha, oldBeta, lap, Q)
        
            [M,N] = size(X); alpha_t = oldAlpha; beta_t = oldBeta; epsilon = 0.01; time=500; e = 100; t = 1;
            
            % start learning iterations
            while e>epsilon && t<time
                
                % E-step
                for s=1:M
                    sample=X(s,:);
                    [estimatedPhi,estimatedGama] = obj.bce_e(alpha_t, beta_t, sample);
                    phiAll(:,:,s) = estimatedPhi;
                    gamaAll(:,s) = estimatedGama;
                end
                
                % M-step
                [alpha_tt,beta_tt] = obj.bce_m(alpha_t, phiAll, gamaAll, X, Q, lap);
                
                % error
                upvalue=0;downvalue=0;
                for index=1:length(Q)
                    upvalue=upvalue+sum(sum(abs(beta_t{index}-beta_tt{index})));
                    downvalue=downvalue+sum(sum(beta_t{index}));
                end
                e=upvalue/downvalue;
%                 disp(['t=',int2str(t),', error=',num2str(e)]);
                
                % update
                alpha_t=alpha_tt;
                beta_t=beta_tt;
                
                t=t+1;
                
            end
            
            resultAlpha=alpha_t;
            resultBeta=beta_t;
            
        end
        
        
        function [phi_t,gama_t] = bce_e(obj, alpha, beta, x)
            
            k = length(alpha); N = size(x,2); V = length(find(x~=0)); filter = ones(k,1)*(x~=0);
            phi_t = ones(k,N)/k.*filter; gama_t = alpha+V/k;
            epsilon = 0.01; time = 500; e = 100; t = 1;
            
            for i=1:k
                for n=1:N
                    if x(n)~=0
                        tempBeta(i,n) = beta{n}(i,x(n));
                    else
                        tempBeta(i,n)=-1;
                    end
                end
            end
            
            while e>epsilon && t<time
                % new phi
                phi_tt=exp((psi(gama_t)-psi(sum(gama_t)))*ones(1,N)).*tempBeta;
                phi_tt=phi_tt./(ones(k,1)*sum(phi_tt+realmin,1));
                phi_tt=phi_tt.*filter;
                
                % new gamma
                gama_tt=alpha+sum(phi_tt,2);
                
                % error of the iteration
                e1=sum(sum(abs(phi_tt-phi_t)))/sum(sum(phi_t));
                e2=sum(abs(gama_tt-gama_t))/sum(gama_t);
                e=max(e1,e2);
                
                % update the variational parameters
                phi_t=phi_tt;
                gama_t=gama_tt;
                % disp(['t=',int2str(t),', e1,e2,e:',num2str(e1),',',num2str(e2),',',num2str(e)]);
                t=t+1;
            end
            
            
        end
        
        
        function [alpha,beta] = bce_m(obj, alpha, phi, gama, X, Q, lap)
            [k,N,M]=size(phi);
            
            for ind=1:N
                beta{ind}=zeros(k,Q(ind));
            end
            
            
            for ind=1:N
                for q=1:Q(ind)
                    temp=zeros(k,N);
                    for s=1:M
                        x=X(s,:);
                        filter=(ones(k,1)*(x==q));
                        temp=temp+phi(:,:,s).*filter;
                    end
                    beta{ind}(:,q)=temp(:,ind);
                end
            end
            
            % smoothing
            for ind=1:N
                beta{ind}=beta{ind}+lap;
                beta{ind}=beta{ind}./(sum(beta{ind},2)*(ones(1,Q(ind)))) ;
            end
            
            alpha_t=alpha;
            epsilon=0.001;
            time=500;
            
            t=0;
            e=100;
            psiGama=psi(gama);
            psiSumGama=psi(sum(gama,1));
            while e>epsilon&&t<time
                g=sum((psiGama-ones(k,1)*psiSumGama),2)+M*(psi(sum(alpha_t))-psi(alpha_t));
                h=-M*psi(1,alpha_t);
                z=M*psi(1,sum(alpha_t));
                c=sum(g./h)/(1/z+sum(1./h));
                delta=(g-c)./h;
                
                % line search
                eta=1;
                alpha_tt=alpha_t-delta;
                while (~isempty(find(alpha_tt<=0, 1)))
                    eta=eta/2;
                    alpha_tt=alpha_t-eta*delta;
                end
                e=sum(abs(alpha_tt-alpha_t))/sum(alpha_t);
                
                alpha_t=alpha_tt;
                
                t=t+1;
            end
            alpha=alpha_t;
        end
        
    end
end