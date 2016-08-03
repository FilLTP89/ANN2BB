function [varargout] = trann_define_inout(varargin)
    
    TnC = varargin{1};
    
    switch TnC
        case 0.5
            % ANN with corner period = 0.5 s
            inp.vTn = [0.6:0.1:1.0,1.25:0.25:5.0];
            tar.vTn = [0,0.05,0.1:0.1:0.5];
        case 0.6
            % ANN with corner period = 0.6 s
            inp.vTn = [0.7:0.1:1.0,1.25:0.25:5.0];
            tar.vTn = [0,0.05,0.1:0.1:0.6];
            
        case 0.75
            % ANN with corner period = 0.75 s
            inp.vTn = [0.8:0.1:1.0,1.25:0.25:5.0];
            tar.vTn = [0,0.05,0.1:0.1:0.7,0.75];
            
        case 1.0
            % ANN with corner period = 1.0 s
            inp.vTn = [1.25:0.25:5.0];
            tar.vTn = [0,0.05,0.1:0.1:1.0];
    end
    
    inp.nT  = length(inp.vTn);
    tar.nT = length(tar.vTn);
    
    varargout{1} = inp.vTn(:);
    varargout{2} = tar.vTn(:);
    varargout{3} = inp.nT;
    varargout{4} = tar.nT;

end