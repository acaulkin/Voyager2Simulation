%Andrew Caulkins
%AERO 343
%Course Project Part 3
%Classical Orbital Elements to a State Vector
%__________________________________________________________________________


%The purpose of this code is to be able to take a vector of the Classical
%Orbital Elements of a satellite (These elements include: 
%semi-major axis [a], eccentricity[e],inclination [inc],
%right ascension of the ascending node [RAAN],
%Argument of Periapsis [AoP], and true anomaly [theta])
%and then perform the necesarry calculations on said satellite
%in order to output the State Vector associated with said Classical Orbital
%Elements

function X = COE2RV(COE, mu)



%Calculate orbital angular momentum from a,mu, and e:
h = sqrt(mu*COE(1)*(1 - COE(2)^2));

%Obtain the Position and Velocity vector(s) in the Perifocal Coordinates:
o = (h^2/mu)*(1/(1 + COE(2)*cos(COE(end))))*[cos(COE(end)),sin(COE(end)),0];
o_dot = (mu/h)*[-sin(COE(end)),COE(2) + cos(COE(end)),0];
o = o';
o_dot = o_dot';
%Defining Transformation Matrices:
O_1 = [cos(COE(5)) sin(COE(5)) 0;
       -sin(COE(5)) cos(COE(5)) 0;
       0 0 1];
O_2 = [1 0 0;
       0 cos(COE(3)) sin(COE(3));
       0 -sin(COE(3)) cos(COE(3))];
O_3 = [cos(COE(4)) sin(COE(4)) 0;
       -sin(COE(4)) cos(COE(4)) 0;
       0 0 1];
    
    
%Using Transformation Matrices to calculate matrix of transformation from
%perifocal to geocentric:
Q = O_1*O_2*O_3;
Q = Q^-1;
%Using matrix to calculate final vectors:
r = Q*o;
r_dot = Q*o_dot;


X = [r;r_dot];


end 

