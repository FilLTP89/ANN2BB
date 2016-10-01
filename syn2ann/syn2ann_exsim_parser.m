%% *Parse outcomes of numerical analysis*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _syn2ann_exsim_parser_: function to parse synthetics from exsim
%% INPUT:
% * _mon (monitor structure)_
% * _mtd (metadata structure)_
%% OUTPUT:
% * _exs (structure of numerical simulations)_
function [varargout] = syn2ann_exsim_parser(varargin)
    %% *SET-UP*
    exs.mon = varargin{1};
    exs.mtd = varargin{2};
    
    for i_ = 1:exs.mon.na
        for j_ = 1:exs.mon.nc
            % file name
            disp(exs.mtd.pf{i_});
            cpp = exs.mon.cp{j_};
            disp(cpp);
            
            str = strcat(exs.mtd.pf{i_},cpp,'_acc_s001.1');
            fid = fopen(str,'r');
            str = textscan(fid,'%f%f','Headerlines',12);
            fclose(fid);
            
            exs.mon.dtm(i_) = mean(diff(str{1}));
            exs.syn{i_}.tha.(cpp) = str{2}./100;
            [~,exs.syn{i_}.thv.(cpp),exs.syn{i_}.thd.(cpp)] = ...
                idc_tha(exs.mon.dtm(i_),exs.syn{i_}.tha.(cpp));
        end
        % time-step number
        exs.mon.ntm(i_) = numel(exs.syn{i_}.tha.(cpp));
        % time-vector
        exs.mon.vtm{i_} = exs.mon.dtm(i_)*(0:exs.mon.ntm(i_)-1);
        
    end
    %% *OUTPUT*
    varargout{1} = exs;
    return
end
