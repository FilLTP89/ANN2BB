function [varargout] = adjust_errbar_tick(varargin)
    erb = varargin{1};
    
    temp = 4:3:length(erb.Xdata);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1;
    % Increase line length by 0.2 units
    erb.Xdata(xleft)  = erb.Xdata(xleft) - .1;
    erb.Xdata(xright) = erb.Xdata(xright) + .1;
    
    varargout{1} = erb;
    return
end