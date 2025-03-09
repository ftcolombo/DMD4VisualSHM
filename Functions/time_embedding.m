function [X1, X2, t] = time_embedding(x, na, t)
% function time_embedding augments the raw data x with na time embeddings

% INPUTS:
% x:        rows are assumed to be measurements, columns are assumed to be time points
% na:       desired number of time embeddings
% t:        time data 

% OUTPUTS:
% X1 & X2:  augmented data matrices, where X2 is X1 shifted one step into the future
% t:        time data 

% F. T. Colombo, September 2024

if na > 1
    Xaug = [];
    for k = 1:na
        Xaug = [Xaug; x(:, k:end-na+k)];
    end
    X1 = Xaug(:, 1:end-1);
    X2 = Xaug(:, 2:end);
else
    X1 = x(:, 1:end-1);
    X2 = x(:, 2:end);
end
t = t(:,1:end-na);
end