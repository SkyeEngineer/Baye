function [X, partition, partition_bin, eta] = irmrnd(alpha, sigma, beta, n)

%
% IRMRND samples from the infinite relational model
%   [X, partition, partition_bin, eta] = irmrnd(alpha, sigma, beta, n)
%
%--------------------------------------------------------------------------
% INPUTS
%   - alpha: scale parameter of the CRP
%   - sigma: discount parameter of the CRP
%   - beta: parameter of the beta distribution for the between-cluster link
%           probabilities
%   - n: Number of nodes
%
% OUPUTS
%   - X: Adjacency matrix (directed case)
%   - partition: partition of the nodes
%   - partition_bin: binary representation of the partition of the nodes
%   - eta: between-clusters link probabilities
%--------------------------------------------------------------------------
% EXAMPLE
% alpha = 1; sigma=0; beta = .1; n=100;
% [X, partition, partition_bin, eta] = irmrnd(alpha, sigma, beta, n)
%--------------------------------------------------------------------------

%
% Copyright Francois Caron, University of Oxford, 2014
%--------------------------------------------------------------------------

% Generate the partition from the (two-parameter) CRP
[partition, ~, K, partition_bin] = pycrprnd(alpha, sigma, n);
% Generate the link probabilities between groups
eta = betarnd(beta, beta, K, K);
% Generate the adjacency matrix
X = rand(n) < partition_bin*eta*partition_bin';



