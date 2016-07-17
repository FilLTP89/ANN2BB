% the subroutine norm_vec.m calculates the  normal fault vector for 
% a given fault geometry 
%
% INPUT PARAMETERS:
% - phi = STRIKE   ( measured clockwise from north.)
% - lambda = RAKE  ( measured counter clockwise from horiz. strike direc.,
%                    as in Aki&Richard )
% - delta =  DIP   ( measured down from the horizontal) 
%
% INPUT PARAMETERS:
% - n = normal fault vector (dimension = 3x1)
%   IMPORTANT:
%       the vecor is given with respet to 
%       the UTM geographic coordinate system (X,Y,Z):
%                X = EAST 
%                Y = NORTH 
%                Z = UP 

function [n] = norm_vec(phi,lambda,delta)

phi = phi*pi/180; % STRIKE (DEG)
lambda = lambda*pi/180; % RAKE (DEG)
delta = delta*pi/180; % DIP (DEG)

% NORMAL FAULT VECTOR n
n = zeros(3,1); 
n(1) = +sin(delta)*cos(phi);
n(2) = -sin(delta)*sin(phi);
n(3) = cos(delta);

return 

  