%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _train_tv_sets_: function to select training and validation percentages
%% *N.B.*
% Need for:_randperm.m_
%% *REFERENCES*
% https://fr.mathworks.com/help/nnet/ug/improve-neural-network-generalization-and-avoid-overfitting.html
% Here a dataset is loaded and divided into two parts: 90% for designing networks and 10% for testing them all.
%
% [x,t] = house_dataset;
% Q = size(x,2);
% Q1 = floor(Q*0.90);
% Q2 = Q-Q1;
% ind = randperm(Q);
% ind1 = ind(1:Q1);
% ind2 = ind(Q1+(1:Q2));
% x1 = x(:,ind1);
% t1 = t(:,ind1);
% x2 = x(:,ind2);
% t2 = t(:,ind2);


function [varargout] = trann_tv_sets(varargin)
    %% *SET-UP*
    nr   = varargin{1};
    pt   = varargin{2};
    
    %% *DEFINE PERCENTAGES*
    Q1   = floor(nr*pt);
    Q2   = nr-Q1;
    idx.all(:,1) = randperm(nr);
    idx.valid    = idx.all(1:Q1,1);
    idx.train    = idx.all(Q1+(1:Q2),1);
    
    %% *OUTPUT*
    varargout{1} = idx.train;
    varargout{2} = idx.valid;
    
    return
end
