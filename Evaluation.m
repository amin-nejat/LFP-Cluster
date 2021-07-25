classdef Evaluation < handle
    properties
        X
    end
    
    methods (Access = private)
      function obj = Evaluation
      end
    end
    
    methods(Static)
      function r = get_instance
         persistent instance;
         if isempty(instance) || ~isvalid(instance)
            instance = Evaluation(); 
         end
         r = instance;
      end
    end
    
    methods
        function [rand_index, f_measure, precision, recall] = clusterers_similarity(obj, A, B)
           if length(A) ~= length(B)
              error('Lengths must be the same');
           end
           
           TP = 0; FP = 0; TN = 0; FN = 0; 
           
           for i = 1 : length(A)
               for j = i+1 : length(A)
                   if A(i) == A(j) && B(i) == B(j)
                       TP = TP + 1;
                   end
                   
                   if A(i) == A(j) && B(i) ~= B(j)
                        FN = FN + 1;
                   end
                   
                   if A(i) ~= A(j) && B(i) == B(j)
                        FP = FP + 1;
                   end
                   
                   if A(i) ~= A(j) && B(i) ~= B(j)
                       TN = TN + 1;
                   end
               end
           end
           rand_index = (TP + TN)/(TP + TN + FP + FN);
           precision = TP / (TP + FP);
           recall = TP / (TP + FN);
           f_measure = 2*precision*recall / (precision + recall);
        end
        
    end
end