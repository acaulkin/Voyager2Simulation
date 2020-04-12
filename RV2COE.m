%Andrew Caulkins
%AERO 343
%Course Project Part 3
%State Vector to Classical Orbital Elements
%__________________________________________________________________________


%The purpose of this code is to be able to take an initial state vector of
%a satellite, and then perform the necesarry calculations on said satellite
%in order to output the Classical Orbital Elements associated with said
%satellite. These elements include: semi-major axis [a], eccentricity[e],
%inclination [inc],right ascension of the ascending node [RAAN],
%Argument of Periapsis [AoP], and true anomaly [theta].

function COE = RV2COE(X, mu)



r = X(1:3);
r_dot = X(4:6);

%% Determining Keplerian Orbital Elements from State Vector:

%First, Calculate the Orbital Momentum Vector:
h_vec = cross(r,r_dot);

%Next, calculate the Eccentricity Vector and the magnitude of it [e]:
e_vec = (cross(r_dot,h_vec)./mu) - r./norm(r);
e = norm(e_vec);
%Determine the Vector n which points towards the ascending node:
n = cross([0,0,1],h_vec); 

%Determining the true anomaly [theta]:
if dot(r,r_dot) >= 0
    theta = acos((dot(e_vec,r))./(norm(e_vec).*norm(r)));
else
    theta = 2*pi - acos((dot(e_vec,r))./(norm(e_vec).*norm(r)));
end

%Determining Orbit Inclination [i]:
inc = acos(h_vec(3)./norm(h_vec));

%Determining Right angle of the Ascensing Node [RAAN]:
if norm(n) == 0
    RAAN = 0;
elseif n(2) >= 0
    RAAN = acos((n(1)./norm(n)));
else
    RAAN = 2*pi - acos((n(1)./norm(n)));
end

%Determing Argument of Periapsis [AoP]:
if dot(n,e_vec) == 0
    AoP = acos(0);
elseif e_vec(3) >= 0
    AoP = acos((dot(n,e_vec))./(norm(n).*norm(e_vec)));
else
    AoP = 2*pi - acos((dot(n,e_vec))./(norm(n).*norm(e_vec)));
end

%Determining the semi-major axis [a]:
a = 1/((2./norm(r)) - ((norm(r_dot).^2)./mu));

%Returning Vector containing all of the Orbital Elements:
COE = [a,e,inc,RAAN,AoP,theta];

end
