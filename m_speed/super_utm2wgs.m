function [varargout] = super_utm2wgs(varargin)
    % This function converts the vectors of UTM coordinates into LAT/LON vectors.
    % Inputs:
    %    E_UTM       - UTM easting in meters
    %    N_UTM       - UTM northing in meters
    %    ZONE_UTM - UTM longitudinal zone
    % Outputs:xx
    %    LON (WGS84 LONgitude vector) in decimal degrees:  ddd.dddddddd
    %    LAT (WGS84 LATitude vector)  in decimal degrees:  ddd.dddddddd
    %
    % Source: DMA Technical Manual 8358.2, Fairfax, VA
    
    %% Argument checking
    narginchk(3,3);         % 3 arguments are required
    E_UTM = varargin{1};
    N_UTM = varargin{2};
    ZONE_UTM = varargin{3};
    
    nE=length(E_UTM);
    nN=length(N_UTM);
    n3=size(ZONE_UTM,1);
    if (nE~=nN || nE~=n3)
        error('x, y and ZONE_UTM vectors should have the same number or rows');
    end
    c=size(ZONE_UTM,2);
    if (c~=3)
        error('ZONE_UTM should be a vector of strings like "30T"');
    end
    
    LON = -999*ones(nE,1);
    LAT = -999*ones(nE,1);
    %% Computing LAT/LON coordinates for each input
    for i=1:nE
        if (ZONE_UTM(i,end)>'X' || ZONE_UTM(i,end)<'C')
            fprintf('utm2wgs: Warning you cannot use lowercase letters in UTM zone\n');
        end
        if (ZONE_UTM(i,end)>'M')
            hemis='N';    % Northern hemisphere
        else
            hemis='S';    % Southern hemisphere
        end
        
        x=E_UTM(i);
        y=N_UTM(i);
        zone=str2double(ZONE_UTM(i,1:2));
        sa = 6378137.000000;                % semi-major axis of the Earth ellipsoid
        sb = 6356752.314245;                % semi-minor axis of the Earth ellipsoid
        e=(((sa^2)-(sb^2))^0.5)/sb;      % squared second eccentricity
        e2= e^2;
        c=sa^2/sb;
        X = x - 500000;
        if hemis=='S' || hemis=='s'
            Y=y-10000000;
        else
            Y=y;
        end
        
        S=((zone*6)-183);
        lat=Y/(6366197.724*0.9996);
        v=(c/((1+(e2*(cos(lat))^2)))^0.5)*0.9996;
        a=X/v;
        a1=sin(2*lat);
        a2=a1*(cos(lat))^2;
        j2=lat+(a1/2);
        j4=((3*j2)+a2)/4;
        j6=((5*j4)+(a2*(cos(lat))^2))/3;
        alpha=(3/4)*e2;
        beta=(5/3)*alpha^2;
        gamma=(35/27)*alpha^3;
        Bm=0.9996*c*(lat-alpha*j2+beta*j4-gamma*j6);
        b=(Y-Bm)/v;
        Epsi=((e2*a^2)/2)*(cos(lat))^2;
        Eps=a*(1-(Epsi/3));
        nab=(b*(1-Epsi))+lat;
        senoheps=(exp(Eps)-exp(-Eps))/2;
        Delta=atan(senoheps/(cos(nab)));
        TaO=atan(cos(Delta)*tan(nab));
        LON(i,1)=(Delta*(180/pi))+S;
        LAT(i,1)=(lat+(1+e2*(cos(lat)^2)-(3/2)*e2*sin(lat)*...
            cos(lat)*(TaO-lat))*(TaO-lat))*(180/pi);
        
    end
    
    %% *OUTPUT*
    varargout{1} = LON;
    varargout{2} = LAT;
    return
end