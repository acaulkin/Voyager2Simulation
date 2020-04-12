function [V1, V2] = LAMBERTSOLVE(R1, R2, DT, opts0) 
% INPUTS:     NAME        DESCRIPTION                      SIZE    TYPE 
%             R1          Position vector at time 1 [km]   3x1     double 
%             R2          Position vector at time 2 [km]   3x1     double 
%             DT          Time between time 1 and 2 [s]    1x1     double 
%             opts0       Options struct with fields 
%                         defined below. 
%             .cBody      Name of central body orbited     NA      string 
%                         ("Sun","Earth","Jupiter", etc.)  
%             .tBody      Name of target body if needed    NA      string 
%             .stopCond   Code stating the stop condition  1x1     int 
%              
% 
% OUTPUTS:    NAME        DESCRIPTION                      SIZE    TYPE 
%             V1          Velocity vector at time 1 [km/s] 3x1     double 
%             V2          Velocity vector at time 2 [km/s] 3x1     double 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mu = GRAVPARAM(opts0.cBody);
    r1norm = norm(R1); r2norm = norm(R2);
    c12 = cross(R1, R2);
    theta = acos(dot(R1,R2)/(r1norm*r2norm));

    if c12(3) <=0
       theta = 2*pi - theta; 
    end
    
    a = sin(theta)*sqrt(r1norm*r2norm/(1-cos(theta)));
    z = -100;
    while F(z,DT, r1norm, r2norm, mu, a) < 0
        z = z + 0.1;
    end

    tol = 1e-10;
    maxiter = 10000;
    ratio = 1;
    n = 0;
    while (abs(ratio) > tol) && (n <=maxiter)
        n = n + 1;
        ratio = F(z,DT, r1norm, r2norm, mu, a)/dFdz(z, r1norm, r2norm,a);
        z = z - ratio;
    end

    f = 1 - y(z, r1norm, r2norm, a)/r1norm;
    g = a*sqrt(y(z, r1norm, r2norm, a)/mu);
    gdot = 1 - y(z, r1norm, r2norm, a)/r2norm;
     
    V1 = 1/g.*(R2 - f.*R1);
    V2 = 1/g.*(gdot.*R2 - R1);
end
function an = y(z, r1norm, r2norm, a)
    an = r1norm + r2norm + a*(z*stumpS(z)-1)/sqrt(stumpC(z));
end
function an = F(z,t, r1norm, r2norm, mu, a)
    an = (y(z, r1norm, r2norm, a)/stumpC(z))^(1.5)*stumpS(z) + a*sqrt(y(z, r1norm, r2norm, a))-sqrt(mu)*t;
end
function an = dFdz(z, r1norm, r2norm,a)
    if z == 0
        an = sqrt(2)/40*y(0, r1norm, r2norm)^1.5 + a/8*(sqrt(y(0, r1norm, r2norm)) + a*sqrt(1/2/y(0, r1norm, r2norm)));
    else
        an = (y(z, r1norm, r2norm, a)/stumpC(z))^1.5*(1/2/z*(stumpC(z) - 3*stumpS(z)/2/stumpC(z))...
            + 3*stumpS(z)^2/4/stumpC(z)) + a/8*(3*stumpS(z)/stumpC(z)*sqrt(y(z, r1norm, r2norm, a))...
            + a*sqrt(stumpC(z)/y(z, r1norm, r2norm, a)));
    end
end
function c = stumpC(z)
%This function evaluates the Stumpff function C(z)
    if z > 0
        c = (1 - cos(sqrt(z)))/z;
    elseif z < 0
        c = (cosh(sqrt(-z)) - 1)/(-z);
    else 
        c = 1/2;
    end
end
function s = stumpS(z)
%This function evaluates the Stumpff function S(z)
    if z > 0
       s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
    elseif z < 0
        s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
    else 
        s = 1/6;
    end
end