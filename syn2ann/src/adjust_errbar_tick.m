function [varargout] = adjust_errbar_tick(varargin)
    erb = varargin{1};
    keyboard
    temp = 4:3:length(erb.YData);
    temp(3:3:end) = [];
    % xleft and xright contain the indices of the left and right
    %  endpoints of the horizontal lines
    xleft = temp; xright = temp+1;
    % Increase line length by 0.2 units
    erb.YData(xleft)  = erb.YData(xleft) - 10;
    erb.YData(xright) = erb.YData(xright) + 10;
    
    varargout{1} = erb;
    return
end