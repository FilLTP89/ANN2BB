function [i,j] = rotrasl3D(str,dip,L,W)
    
    %=================== function rotrasl3D ============================%
    %
    % Given the geometry of the extended fault
    %         INPUT PARAMETERS
    % strike: str (in deg!)
    % dip: dip (in deg!)
    % length (along strike): L (m/km: as you like!)
    % width (along dip): W (m/km: as you like!)
    %
    % the function computes:
    % 1. 3D rotation of the fault plan 1-2-3-4 (see below):
    %       1a. about the vertical axis 1-4 of angle=str
    %       1b. about the axis through points 1-2 of angle = dip
    %
    % --------------------------------
    %       |  Reference frame  |
    % --------------------------------
    %             ^ Z (Up)
    %             |
    %             |
    %  origin <-- 1------------2 ---> X (East)
    %             |     P      |
    %             |            |
    %             4------------3
    % 1 coincides with 0 = origin of the ref. system.
    %
    % OUTPUT:
    % - versors v1 and v2 (size: 3x1) which allow the transformation of any
    % point P on the fault plane according to the performed rotation
    % (strike-dip).
    %
    % NOTE1: the following geographic reference system (UTM: X = EAST, Y = NORTH)
    % is adopted:
    %
    %   o = pointing up
    %   x = pointing down
    %
    %   -----------------------
    %   |  Reference system |
    %   -----------------------
    %   Z = UP
    %   ^
    %   |
    %   |
    %   |
    %   | Y (NS)
    %   o------------>  X  = EW
    %
    % NOTE2: uses the function roto3D.m
    %===================================================================%
    
    % Defines corners of the fault
    %   1------------2
    %   |            |
    %   |            |
    %   4------------3
    % length 12 = L
    % length 14 = W
    
    p = zeros(3,4);
    p(:,1) = [0 0 0]';
    p(:,2) = [L 0 0]';
    p(:,3) = [L 0 W]';
    p(:,4) = [0 0 W]';
    
    x = zeros(4,1);
    y = zeros(4,1);
    z = zeros(4,1);
    for j = 1:5
        x(:) = p(1,:);
        y(:) = p(2,:);
        z(:) = p(3,:);
    end
    
    % 1a) rotation according to str (vertical axis z through points 1,4)
    % effective angle of rotation phi_e = -(str - pi/2)
    % uses function roto3D.m : [vp] = roto3D(v,axis,theta);
    ax = 3;
    str = str*pi/180;
    phi_e = - (str - pi/2);
    p2 = zeros(3,4);
    for j=1:4
        [p2(:,j)] = roto3D(p(:,j),ax,phi_e);
    end
    % updates coordinates x,y,z
    x(:) = p2(1,:);
    y(:) = p2(2,:);
    z(:) = p2(3,:);
    
    % 1b) rotation according to dip (horizontal axis through points 1,2
    %     i.e. along strike)
    %
    % Define orthonormal system along rotation axis (along strike direction)
    vx=0;
    vy=0;
    m=abs(z(3)-z(2));
    vz=(z(3)-z(2))/m;
    %
    ix=x(2)-x(1);
    iy=y(2)-y(1);
    iz=0;
    n=sqrt(ix*ix+iy*iy);
    ix=ix/n;
    iy=iy/n;
    %
    wx=iy*vz-iz*vy;
    wy=iz*vx-ix*vz;
    wz=ix*vy-iy*vx;
    %
    dip = -dip*pi/180;
    ux=(cos(dip)*wx+sin(dip)*vx)*m;
    uy=(cos(dip)*wy+sin(dip)*vy)*m;
    uz=(cos(dip)*wz+sin(dip)*vz)*m;
    x(3)=x(2)+ux;
    y(3)=y(2)+uy;
    z(3)=z(2)+uz;
    x(4)=x(1)+ux;
    y(4)=y(1)+uy;
    z(4)=z(1)+uz;
    
    
    i(1)=x(2)-x(1);
    i(2)=y(2)-y(1);
    i(3)=z(2)-z(1);
    i=i/norm(i);
    j(1)=x(4)-x(1);
    j(2)=y(4)-y(1);
    j(3)=z(4)-z(1);
    j=j/norm(j);
    
    return
end