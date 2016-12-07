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
%
    % slip_vec.m calculates the slip vector for a given fault geometry
    %
    % INPUTS:
    % - phi = STRIKE   ( measured clockwise from north.)
    % - lambda = RAKE  ( measured counter clockwise from horiz. strike direc.,
    %                    as in Aki&Richard )
    % - delta =  DIP   ( measured down from the horizontal)
    % - varargin = reference system convention
    % OUTPUTS:
    % - slp = slip vector (dimension = 3x1)
    %   N.B.:
    %   if ref_sys = 0
    %       the vecor is given with respet to
    %       the UTM geographic coordinate system (X,Y,Z):
    %                X = EAST
    %                Y = NORTH
    %                Z = UP
    %   if ref_sys = 1
    %       the vecor is given with respet to
    %       the HISADA geographic coordinate system (X,Y,Z):
    %                X = NORTH
    %                Y = EAST
    %                Z = DOWN
function [varargout] = slip_vec(varargin)
    %% *SET-UP*
    phi    = varargin{1}*pi/180; % STRIKE (DEG)
    lambda = varargin{2}*pi/180; % RAKE (DEG)
    delta  = varargin{3}*pi/180; % DIP (DEG)
    ref_sys = 0;
    if nargin>3
        ref_sys = varargin{4};
    end
    
    %% *SLIP VECTOR slp*
    slp = zeros(3,1);
    if ref_sys==0
        disp('UTM REFSYS');
        slp(1) = +cos(lambda)*sin(phi) - sin(lambda)*cos(delta)*cos(phi);
        slp(2) = +cos(lambda)*cos(phi) + sin(lambda)*cos(delta)*sin(phi);
        slp(3) = sin(lambda)*sin(delta);
    elseif ref_sys==1
        disp('HISADA REFSYS');
        slp(1) = +cos(lambda)*cos(phi) + cos(delta)*sin(lambda)*sin(phi);
        slp(2) = +cos(lambda)*sin(phi) - cos(delta)*sin(lambda)*cos(phi);
        slp(3) = -sin(lambda)*sin(delta);
    end
    
    %% *OUTPUT*
    varargout{1} = round(slp.*1e8).*1e-8;
    return
end