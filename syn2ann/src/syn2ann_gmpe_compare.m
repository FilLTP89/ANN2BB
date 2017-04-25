%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_gmpe_compare_: function to compare SM indicators with GMPEs
%%
%% *N.B.*
function [varargout]=syn2ann_gmpe_compare(varargin)
    
    % default inputs
    def.mod   = 'lanzano_et_al_2016_rfocal';
    def.Mw   = -999;
    def.MJMA = -999;
    def.Delta = -999;
    def.hrupt = -999;
    def.lons  = -999;
    def.lats  = -999;
    def.fsof  = 'UN';
    def.EC8   = 'A';
    def.bas   = 0;
    def.fld   = ['/home/filippo/Data/Filippo/PHD_passing_through_polimi/',...
        'syn2ann/database/maps'];
    def.map   = {'null'};
    def.var   = {'pga'};
    def.cpn   = {'Rgh [m/s2]'};
    def.hyp   = zeros(3,1);
    def.utmzone = 'T32';
    
    % parser parameters
    inp = inputParser;
    addParameter(inp,'mod'  ,def.mod  ,@ischar);
    addParameter(inp,'Mw'  ,def.Mw  ,@isnumeric);
    addParameter(inp,'MJMA',def.MJMA,@isnumeric);
    addParameter(inp,'hrupt',def.hrupt,@isnumeric);
    addParameter(inp,'fsof',def.fsof,@ischar);
    addParameter(inp,'EC8',def.EC8,@ischar);
    addParameter(inp,'bas',def.bas,@isnumeric);
    addParameter(inp,'map',def.map,@iscell);
    addParameter(inp,'var',def.var,@iscell);
    addParameter(inp,'fld',def.fld,@ischar);
    addParameter(inp,'cpn',def.cpn,@iscell);
    addParameter(inp,'hyp',def.hyp,@isnumeric);
    addParameter(inp,'utmzone',def.utmzone,@ischar);
    
    % parse results
    parse(inp,varargin{:});
    
    
    for i_=1:numel(inp.Results.map)
        map = importdata(fullfile(inp.Results.fld,inp.Results.map{i_}));
        [~,idx.ew,~] = intersect(map.textdata,'E-UTM [m]');
        [~,idx.ns,~] = intersect(map.textdata,'N-UTM [m]');
        [~,idx.ud,~] = intersect(map.textdata,'U-UTM [m]');
        crd.ew = map.data(:,idx.ew);
        crd.ns = map.data(:,idx.ns);
        if isempty(idx.ud)
            crd.ud = zeros(size(crd.ew));
        else
            crd.ud = map.data(:,idx.ud);
        end
        for j_=1:numel(inp.Results.var)
            [~,idx.var(j_),~] = intersect(map.textdata,inp.Results.cpn{j_});
        end
        % _focal distance_
        
        Delta = [crd.ew,crd.ns,crd.ud]-...
            repmat(inp.Results.hyp(:)',[numel(crd.ew),1]);
        Delta = sqrt(sum(Delta.*Delta,2))./1e3;
        [crd.lon,crd.lat] = super_utm2wgs(crd.ew,crd.ns,...
            repmat(inp.Results.utmzone,[numel(crd.ew),1]));
        
        if numel(inp.Results.Mw)==1
            Mw=repmat(inp.Results.Mw,[numel(Delta),1]);
        else
            Mw=inp.Results.Mw;
        end
        if numel(inp.Results.EC8)==1
            EC8 = repmat(inp.Results.EC8,[numel(Delta),1]);
        else
            EC8 = inp.Results.EC8;
        end
        if numel(inp.Results.bas)==1
            bas = repmat(inp.Results.bas,[numel(Delta),1]);
        else
            bas = inp.Results.bas;
        end
        
        %% *GMPE*
        gmpe = generate_gmpe('mod',inp.Results.mod,...
            'Mw',Mw,...
            'Delta',Delta,...
            'hrupt',inp.Results.hrupt,...
            'lons',crd.lon,...
            'lats',crd.lat,...
            'fsof',inp.Results.fsof,...
            'EC8',EC8,'bas',bas);
        
        keyboard
    end
    
    %% *SIMULATIONS*
    return
end