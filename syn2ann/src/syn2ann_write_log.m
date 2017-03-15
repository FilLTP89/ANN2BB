function syn2ann_write_log(sas,varargin)
    %% *SET-UP*
    def.pfn = './prova.csv';
    def.stk = 0;
    % parser parameters
    inp = inputParser;
    addParameter(inp,'pfn',def.pfn,@ischar);
    addParameter(inp,'stk',def.stk,@isnumeric);
    
    % parse results
    parse(inp,varargin{:});
    pfn = inp.Results.pfn;
    
    %% *WRITE LOG FILE*
    for i_=1:sas.mon.na
        E{i_,1} = sas.mon.eutm(i_);
        N{i_,1} = sas.mon.nutm(i_);
        for j_=1:sas.mon.nc
            pga.(sas.mon.cp{j_}){i_,1} = abs(sas.syn{i_}.pgv.(sas.mon.cp{j_})(2));
        end
        [pga.pl{i_,1},pga.pp{i_,1}] = rotationComponent2D(sas.syn{i_}.pgv.ew(2),...
            sas.syn{i_}.pgv.ns(2),inp.Results.stk);
        pga.pl{i_,1} = abs(pga.pl{i_,1});
        pga.pp{i_,1} = abs(pga.pp{i_,1});
        pga.gh{i_,1} = sqrt(abs(sas.syn{i_}.pgv.ew(2)*sas.syn{i_}.pgv.ns(2)));
    end
    super_csvwrite(inp.Results.pfn,...
        [E,N,pga.ew,pga.ns,pga.gh,pga.pl,pga.pp],...
        '%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f\n',...
        {'E-UTM - m','N-UTM - m','PGV-EW - m/s2','PGV-NS - m/s2','PGV-GH - m/s2','PGV-FP - m/s2','PGV-FO - m/s2'});
            
    
    return
end
