%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _test_ann_justPSA_: function train ANN on PSA values
%% *N.B.*
% Need for:
% _trann_define_inout.m, trann_check_vTn.m, trann_tv_sets.m_
% _ANN MATLAB tool_
%% *REFERENCES*
function [varargout] = test_ann_justPSA(varargin)
    %% *SET-UP*
    ann = varargin{1};
    rec = varargin{2};
    
    inn = cell(rec.mon.nc,1);
    %
    % _define input/target natural periods_
    %
    [inp.vTn,tar.vTn,inp.nT,tar.nT] = trann_define_inout(ann.TnC);
    [inp.idx,tar.idx] = trann_check_vTn(inp,tar,rec.mon,1e-8);
    ann.mon.vTn = [tar.vTn(:);inp.vTn(:)];
    
    for i_=1:rec.mon.na
        %
        % _check input/target natural periods with database_
        %
        flag.ann = seismo_dir_conversion(ann.cpp);
        inp.psa = ones(numel(inp.idx),1);
        %         if strcmpi(ann.cpp,'gh')
        fprintf('TESTING ANN ON COMPONENT:\n');
        count=0;
        for j_ = 1:rec.mon.nc
            cpp.rec = rec.mon.cp{j_};
            flag.rec = seismo_dir_conversion(cpp.rec);
            
            if any(strcmpi(flag.rec,flag.ann))
                fprintf('---> %s FOUND!\n',cpp.rec);
                count = count+1;
                inp.psa = inp.psa.*rec.syn{i_}.psa.(cpp.rec)(inp.idx,1);
            end
        end
        inp.psa = inp.psa.^(1/count);
        out.psa = 10.^(sim(ann.net,log10(inp.psa*100)));
        ann.syn{i_}.psa.(ann.cpp) = [out.psa(:)./100;inp.psa(:)];
    end
    
    %% *OUTPUT*
    varargout{1} = ann;
    
    return
end
