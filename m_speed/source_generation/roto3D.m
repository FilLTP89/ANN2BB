%**************************************************************************
% Copyright (C) 2012 The SPEED FOUNDATION
% Author: SPEED-Team (POLIMI)
%         Politecnico di Milano 
%         P.zza Leonardo da Vinci, 32 
%         20133 Milano 
%         Italy                                        
%
% This file is part of SPEED.
%
% SPEED is free software; you can redistribute it and/or modify it
% under the terms of the GNU Affero General Public License as
% published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
%
% SPEED is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with SPEED.  If not, see <http://www.gnu.org/licenses/>.
%**************************************************************************



function [vp] = roto3D(v,axis,theta)
% 3D ROTATION OF ANGLE THETA about either x,y or z direction 
%
% INPUTS
% v = 3-component vector to be rotated of angle theta 
%     IMPORTANT: v must have dimension: (3x1)!
% axis = axis of rotation (1=x;2=y;3=z)
% theta = angle of rotation (in rad) 
%     IMPORTANT: theta measured positive in the counterclockwise direction!
%
% OUTPUTS 
% vp = rotated vector of dimension (3x1)
%
% REFERENCE: Foley & Van Dam, Chapter 5 (Lecture_07_6.pdf )


R = zeros(4,4); 
v = [v;1]; 

if axis == 1
    R = [1    0          0       0;
        0 cos(theta) -sin(theta) 0;
        0 sin(theta) cos(theta)  0;
        0     0           0      1];
elseif axis == 2
    R = [cos(theta) 0 -sin(theta) 0;
        0           1    0        0;
        sin(theta)  0 cos(theta)  0;
        0     0           0      1];
elseif axis == 3
    R = [cos(theta) -sin(theta) 0 0;
        sin(theta) cos(theta)   0 0;
        0     0          1        0;
        0     0          0       1];    
end

vp = R*v; 
vp = vp(1:3); 
return 
