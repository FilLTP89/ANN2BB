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

function [Xp] = fault_points(v1,v2,X)

% Output: x-y points on the fault plane

Xp = zeros(2,1); 

P = [X(1) X(2) X(3)]; 
Xp(1) = (v2(1)*P(2)-v2(2)*P(1))/(v1(2)*v2(1)-v2(2)*v1(1));
Xp(2) = -(P(3)-v1(3)*Xp(1))/v2(3); 

return

