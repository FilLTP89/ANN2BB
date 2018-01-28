function [amax] = newmark_sa(gacc,T,zeta,dt)
    %[amax] = newmark_sa(gacc,T,zeta,dt)
    %   Compute response spectrum via Newmark method for direct integration
    %   INPUT:
    %       * gacc: ground acceleration [npun x 1]
    %       * T: structural period (s)
    %       * zeta: damping ratio
    %       * dt: time step for numerical integration
    %   OUTPUT:
    %       * ymax: maximum relative displacement of single dof
    
    % newmark coefficients
    beta = 0.25;
    gamma = 0.5;
    
    % initial conditions
    y0 = 0;
    yp0 = 0;
    
    % initialization
    npun = length(gacc);
    y = zeros(npun,1);
    yp = zeros(npun,1);
    ypp = zeros(npun,1);
    
    % natural circular frequency of sdof system
    wn = 2*pi/T;
    
    
    y(1) = y0;
    yp(1) = yp0;
    ypp(1) = -gacc(1)-2*wn*zeta*yp0-wn^2*y0;
    
    % Integration coefficients
    keff = wn^2 + 1/(beta*dt^2) + gamma*2*wn*zeta/(beta*dt);
    a1 = 1/(beta*dt^2)+gamma*2*wn*zeta/(beta*dt);
    a2 = 1/(beta*dt)+2*wn*zeta*(gamma/beta-1);
    a3 = (1/(2*beta)-1)+2*wn*zeta*dt*(gamma/(2*beta)-1);
    
    
    for i=1:npun-1
        y(i+1)   = (-gacc(i+1)+a1*y(i)+a2*yp(i)+a3*ypp(i))/keff;
        ypp(i+1) = (y(i+1)-y(i)-dt*yp(i)-dt^2*ypp(i)/2)/(beta*dt^2) + ypp(i);
        yp(i+1)  = yp(i)+dt*ypp(i)+dt*gamma*(ypp(i+1)-ypp(i));
    end
    
    ymax = max(abs(y));
    amax = ymax*(2*pi/T)^2;
    
    return
end
