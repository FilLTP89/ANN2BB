function [varargout] = trann_define_inout_reduced(varargin)
    
    TnC = varargin{1};
    
    switch TnC
        case 0.5
            % ANN with corner period = 0.50 s
            inp.vTn = [0.5:0.1:1.0,1.25:0.25:2.00];
            tar.vTn = [0,0.05,0.1:0.1:0.4];
        case 0.6
            % ANN with corner period = 0.60 s
            inp.vTn = [0.6:0.1:1.0,1.25:0.25:2.00];
            tar.vTn = [0,0.05,0.1:0.1:0.5];
            
        case 0.75
            % ANN with corner period = 0.75 s
            inp.vTn = [0.75,0.8:0.1:1.0,1.25:0.25:2.00];
            tar.vTn = [0,0.05,0.1:0.1:0.7];
        case 1.0
            % ANN with corner period = 1.00 s
            inp.vTn = [1.0:0.25:2.00];
            tar.vTn = [0,0.05,0.1:0.1:0.9];
        case 0.25
            % ANN with corner period = 0.25 s
            inp.vTn = [0.25:0.25:2.00];
            tar.vTn = [0:0.05:0.20];
            
    end
    
    inp.nT = length(inp.vTn);
    tar.nT = length(tar.vTn);
    
    varargout{1} = inp.vTn(:);
    varargout{2} = tar.vTn(:);
    varargout{3} = inp.nT;
    varargout{4} = tar.nT;
    
end

