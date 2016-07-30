function [varargout] = trann_tv_sets(varargin)

    nr   = varargin{1};
    pt   = varargin{2};
    Q1   = floor(nr*pt);
    Q2   = nr-Q1;
    idx.all(:,1) = randperm(nr);
    idx.valid    = idx.all(1:Q1,1);
    idx.train    = idx.all(Q1+(1:Q2),1);
    %idx.train    = randi(nr,nr,1);
    %idx.valid    = randi(nr,round(pt*nr),1);

    %% *OUTPUT*
    varargout{1} = idx.train;
    varargout{2} = idx.valid;

return
end
