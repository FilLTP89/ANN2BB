%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2amon.na_write_maps_: function to write contour maps
%% *N.B.*
% Need for:
% __

if mon.map.flg
    fprintf('============================\n');
    fprintf('----------1. WRITE MAPS-----\n');
    fprintf('============================\n');
    
    [map(1:mon.na).Geometry] = deal('Point');
    for i_ = 1:mon.na
        [map(i_).X] = mon.eutm(i_);
        [map(i_).Y] = mon.nutm(i_);
    end
    
    %% *WRITE PGA MAPS*
    if any(strcmpi(mon.map.typ,'pga'))
        for i_=1:mon.na
            for j_=1:mon.nc
                map(i_).(strcat('R',mon.cp{j_})) = ...
                    abs(spm.sps.(mon.cp{j_}).syn{i_}.pga.(mon.cp{j_})(2));
            end
            [map(i_).Rp,map(i_).Ro] = ...
                rotationComponent2D(spm.sps.ew.syn{i_}.pga.ew(2),...
                spm.sps.ns.syn{i_}.pga.ns(2),mon.map.stk);
            map(i_).Rp = abs(map(i_).Rp);
            map(i_).Ro = abs(map(i_).Ro);
            map(i_).Rgh = ...
                sqrt(abs(spm.sps.ew.syn{i_}.pga.ew(2)*spm.sps.ns.syn{i_}.pga.ns(2)));
        end
        shapewrite(map,strcat(mon.map.fnm,'_pga'));
    end
    
    %% *WRITE PGV MAPS*
    if any(strcmpi(mon.map.typ,'pgv'))
        for i_=1:mon.na
            for j_=1:mon.nc
                map(i_).(strcat('R',mon.cp{j_})) = ...
                    abs(spm.sps.(mon.cp{j_}).syn{i_}.pgv.(mon.cp{j_})(2));
            end
            [map(i_).Rp,map(i_).Ro] = ...
                rotationComponent2D(spm.sps.ew.syn{i_}.pgv.ew(2),...
                spm.sps.ns.syn{i_}.pgv.ns(2),mon.map.stk);
            map(i_).Rp = abs(map(i_).Rp);
            map(i_).Ro = abs(map(i_).Ro);
            map(i_).Rgh = ...
                sqrt(abs(spm.sps.ew.syn{i_}.pgv.ew(2)*spm.sps.ns.syn{i_}.pgv.ns(2)));
        end
        shapewrite(map,strcat(mon.map.fnm,'_pgv'));
    end
    
    %% *WRITE PGD MAPS*
    if any(strcmpi(mon.map.typ,'pgd'))
        for i_=1:mon.na
            for j_=1:mon.nc
                map(i_).(strcat('R',mon.cp{j_})) = ...
                    abs(spm.sps.(mon.cp{j_}).syn{i_}.pgd.(mon.cp{j_})(2));
            end
            [map(i_).Rp,map(i_).Ro] = ...
                rotationComponent2D(spm.sps.ew.syn{i_}.pgd.ew(2),...
                spm.sps.ns.syn{i_}.pgd.ns(2),mon.map.stk);
            map(i_).Rp = abs(map(i_).Rp);
            map(i_).Ro = abs(map(i_).Ro);
            map(i_).Rgh = ...
                sqrt(abs(spm.sps.ew.syn{i_}.pgd.ew(2)*spm.sps.ns.syn{i_}.pgd.ns(2)));
        end
        shapewrite(map,strcat(mon.map.fnm,'_pgd'));
    end
    
    %% *WRITE RSD MAPS*
    if any(strcmpi(mon.map.typ,'rsd'))
        idx = find(spm.sps.ew.mon.vTn==mon.map.vTn.rsd);
        for i_=1:mon.na
            for j_=1:mon.nc
                map(i_).(strcat('R',mon.cp{j_})) = spm.sps.(mon.cp{j_}).syn{i_}.rsd.(mon.cp{j_})(idx);
            end
            map(i_).Rgh = sqrt(spm.sps.ew.syn{i_}.rsd.ew(idx)*spm.sps.ns.syn{i_}.rsd.ns(idx));
        end
        shapewrite(map,strcat(mon.map.fnm,'_rsd'));
    end
    
    %% *WRITE PSA MAPS*
    if any(strcmpi(mon.map.typ,'rsd'))
        idx = find(spm.sps.ew.mon.vTn==mon.map.vTn.psa);
        for i_=1:mon.na
            for j_=1:mon.nc
                map(i_).(strcat('R',mon.cp{j_})) = ...
                    abs(spm.sps.(mon.cp{j_}).syn{i_}.psa.(mon.cp{j_})(idx));
            end
            map(i_).Rgh =...
                sqrt(abs(spm.sps.ew.syn{i_}.psa.ew(idx)*spm.sps.ns.syn{i_}.psa.ns(idx)));
        end
        shapewrite(map,strcat(mon.map.fnm,'_psa'));
    end
end