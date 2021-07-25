classdef (Sealed) Dimensionality_Reduction < handle
    properties
        X
        new_dim
    end
    
    methods (Access = private)
      function obj = Dimensionality_Reduction
      end
   end

    methods(Static)
      function r = get_instance
         persistent instance;
         if isempty(instance)
            instance = Dimensionality_Reduction; 
         end
         r = instance;
      end
    end
    
    methods
        function reduced = pca_reduction(obj)
            [coeff, score] = pca(obj.X);
            reduced = obj.X * coeff(:,1:obj.new_dim);
            
        end
        
        function reduced = svd_reduction(obj)
            [U, S, V] = svd(obj.X);
            reduced = U(:, 1:k)*S(1:k, 1:k);
        end
        
        function pca_visualize(obj)
            mapcaplot(obj.X);
        end
        
        function svd_getk(obj)
            [U, S, V] = svd(obj.X);
            costs = diag(S);
            figure(1);
            plot(linspace(1, length(costs), length(costs)), costs, 'color', 'blue');
            hold on;
        end
        
        function kmeans_getk(obj)
            costs = zeros(1, 14);
            for k = 2 : 15
                [idx, C, sumd] = kmeans(obj.X, k);
                costs(k-1) = sum(sumd);
            end
            figure(1);
            plot(linspace(2, 15, length(costs)), costs, 'color', 'blue');
            hold on;
        end
    end
end