% the subroutine slip_vec.m calculates the slip vector for a given fault
% geometry 
%
% INPUT PARAMETERS:
% - phi = STRIKE   ( measured clockwise from north.)
% - lambda = RAKE  ( measured counter clockwise from horiz. strike direc.,
%                    as in Aki&Richard )
% - delta =  DIP   ( measured down from the horizontal) 
%
% INPUT PARAMETERS:
% - s = slip vector (dimension = 3x1)
%   IMPORTANT:
%       the vecor is given with respet to 
%       the UTM geographic coordinate system (X,Y,Z):
%                X = EAST 
%                Y = NORTH 
%                Z = UP 

function [s] = slip_vec(phi,lambda,delta); 

phi = phi*pi/180; % STRIKE (DEG)
lambda = lambda*pi/180; % RAKE (DEG)
delta = delta*pi/180; % DIP (DEG)

% SLIP VECTOR s
s = zeros(3,1);
s(1) = +cos(lambda)*sin(phi) - sin(lambda)*cos(delta)*cos(phi);
s(2) = +cos(lambda)*cos(phi) + sin(lambda)*cos(delta)*sin(phi); 
s(3) = sin(lambda)*sin(delta); 

return 

  