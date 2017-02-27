%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_group_map_files_: function to read csv map files and group them
% in a unique csv file (or shape file if the toolbox is installed)
%% *N.B.*
% Need for:
% __
ccc;
maps.pfn = '/media/filippo/Data/Filippo/PHD_passing_through_polimi/syn2ann/database/maps/maps_mrn';
maps.tag.org = 'map_mrn';
maps.tag.all = 'maps_ppe2012';
maps.typ = {'spm';'hyb'};
maps.cnt = {'pga';'pgv';'psa_10';'psa_20';'psa_50'};
maps.fid = 1:100;
maps.nb.typ = numel(maps.typ);
maps.nb.cnt = numel(maps.cnt);
maps.nb.fid = numel(maps.fid);

for i_=1:maps.nb.typ
    for j_=1:maps.nb.cnt
        fnm.all = fullfile(maps.pfn,strcat(maps.tag.all,'_',maps.typ{i_},'_',maps.cnt{j_},'.csv'));
        fid.all = fopen(fnm.all,'w+');
        for k_=1:1
            fnm.sng{j_} = fullfile(maps.pfn,strcat(maps.tag.org,'_',maps.typ{i_},...
                '_',maps.cnt{j_},'_',num2str(maps.fid(k_)),'.csv'));
            maps.hln = importdata(fnm.sng{j_});
            maps.hln = maps.hln.textdata;
            maps.nb.hln = numel(maps.hln);
            for ii_=1:maps.nb.hln
                fprintf(fid.all,'%s,',maps.hln{ii_});
            end
        end
        fprintf(fid.all,'\n');
        for k_=1:maps.nb.fid
            fnm.sng{j_} = fullfile(maps.pfn,strcat(maps.tag.org,'_',maps.typ{i_},...
                '_',maps.cnt{j_},'_',num2str(maps.fid(k_)),'.csv'));
            dat = importdata(fnm.sng{j_});
            dat = dat.data;
            for ii_=1:size(dat,1)
                for jj_=1:maps.nb.hln
                    fprintf(fid.all,'%20.6f,',dat(ii_,jj_));
                end
                fprintf(fid.all,'\n');
            end
        end
    end
end