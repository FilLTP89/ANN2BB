function [Xp] = fault_points(v1,v2,X)

Xp = zeros(2,1); % x-y on fault plane

P = [X(1) X(2) X(3)]; 
Xp(1) = (v2(1)*P(2)-v2(2)*P(1))/(v1(2)*v2(1)-v2(2)*v1(1));
Xp(2) = -(P(3)-v1(3)*Xp(1))/v2(3); 

return

