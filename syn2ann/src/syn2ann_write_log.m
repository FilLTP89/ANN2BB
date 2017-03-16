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
        sts{i_,1} = sas.mon.st{i_};
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
        [sts,E,N,pga.ew,pga.ns,pga.gh,pga.pl,pga.pp],...
        '%s,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f,%20.3f\n',...
        {'STATION','E-UTM - m','N-UTM - m','PGV-EW - m/s','PGV-NS - m/s','PGV-GH - m/s','PGV-FP - m/s','PGV-FO - m/s'});
            
    
    return
end
