%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_pbs_drive_: function to parse PBS records
%% *N.B.*
% Need for:
% _syn2ann_pbs_parser.m,syn2ann_thp.m,syn2ann_spp.m_

fprintf('============================\n');
fprintf('------2. PB SIMULATIONS-----\n');
fprintf('============================\n');

%% *PARSING ORIGINAL NUMERICAL SIMULATIONS*
fprintf('--> Parsing\n');
mon.lfr = 0.05;
mon.hfr = [];
flag_aae=1;

[mon,pbs.org]= syn2ann_pbs_parser(mon,bhr);

if flag_aae==1
    for i=1:pbs.org.mon.na
        for j=1:pbs.org.mon.nc
            pbs.org.syn{i}.thd.(pbs.org.mon.cp{j}) = rec.fil.syn{i}.thd.(pbs.org.mon.cp{j}); 
            pbs.org.syn{i}.thv.(pbs.org.mon.cp{j}) = rec.fil.syn{i}.thv.(pbs.org.mon.cp{j}); 
            pbs.org.syn{i}.tha.(pbs.org.mon.cp{j}) = rec.fil.syn{i}.tha.(pbs.org.mon.cp{j}); 
            %[pbs.org.syn{i}.tha.(pbs.org.mon.cp{j}),...
            %    pbs.org.syn{i}.thv.(pbs.org.mon.cp{j}),... 
            %    pbs.org.syn{i}.thd.(pbs.org.mon.cp{j}),... 
            %    pbs.org.mon.vtm{i},~] = ...
            %    bpf_tha(pbs.org.mon.dtm(i),...
            %    rec.fil.syn{i}.tha.(pbs.org.mon.cp{j})(1:2:end),0.05,[]);
        end
        pbs.org.mon.dtm(i) = rec.fil.mon.dtm(i);
        pbs.org.mon.ntm(i) = rec.fil.mon.ntm(i);
        pbs.org.mon.vtm{i} = rec.fil.mon.vtm{i};
    end
end
%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
pbs.org = syn2ann_thp(pbs.org);

%% *SPECTRA*
fprintf('--> Spectra\n');
[pbs.org] = syn2ann_spp(pbs.org);
