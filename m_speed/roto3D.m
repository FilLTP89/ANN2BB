% 3D ROTATION OF ANGLE THETA about either x,y or z direction 

% input parameters:
% v = 3-component vector to be rotated of angle theta 
% IMPORTANT: v must have dimension: (3x1)!
% axis = axis of rotation (1=x;2=y;3=z)
% theta = angle of rotation (in rad) 
% IMPORTANT: theta measured positive in the counterclockwise direction!
%
% output: 
% vp = rotated vector of dimension (3x1)

% REFERENCE: Foley & Van Dam, Chapter 5 (Lecture_07_6.pdf )

function [vp] = roto3D(v,axis,theta)

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
