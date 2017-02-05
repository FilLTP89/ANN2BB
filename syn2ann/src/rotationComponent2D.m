function [vp,vo] = rotationComponent2D(vx,vy,alpha)
% alpha = azimuth measured clockwise from the North (in rad) 
% vx = east-west component 
% vy = north-south component
if alpha ~= 0 
    vp = vx.*sin(alpha) + vy.*cos(alpha); % -> parallel  
    vo = vx.*cos(alpha) - vy.*sin(alpha); % -> orthogonal   
else
    vp = vx;
    vo = vy;
end
        



