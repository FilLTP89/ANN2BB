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

function [n] = norm_vec(varargin)
    phi    = varargin{1};
    lambda = varargin{2};
    delta  = varargin{3};
    %
    % norm_vec.m calculates the  normal fault vector for
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
    
    ref_sys = 0;
    if nargin>3
        ref_sys = varargin{4};
    end
    
    phi    = phi*pi/180;    % STRIKE (DEG)
    lambda = lambda*pi/180; % RAKE (DEG)
    delta  = delta*pi/180;  % DIP (DEG)
    
    % NORMAL FAULT VECTOR n
    n = zeros(3,1);
    if ref_sys==0
        n(1) = +sin(delta)*cos(phi);
        n(2) = -sin(delta)*sin(phi);
        n(3) = +cos(delta);
    elseif ref_sys==1
        n(1) = -sin(delta)*sin(phi);
        n(2) = +sin(delta)*cos(phi);
        n(3) = -cos(delta);
    end
    n = round(n.*1e8).*1e-8;
    return
end