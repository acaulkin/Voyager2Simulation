function Xsc_s = PLANET2SUN(t,Xsc_p,planet)
%%% This function transforms ICRF states from Planet centered to Sun
%%% centered.
%%%
%%% INPUTS:     NAME        DESCRIPTION                      SIZE    TYPE
%%%             Xsc_p       State vector relative to planet  6x1     double
%%%                         [x,y,z,vx,vy,vz]'
%%%                         [km,km,km,km/s,km/s,km/s]'
%%%             t           Time of state in s since ref     1x1     double
%%%             planet      Name of planet for state vector  NA      string
%%%             
%%%
%%% OUTPUTS:    NAME        DESCRIPTION                      SIZE    TYPE
%%%             Xsc_s       State vector relative to Sun     6x1     double
%%%                         [x,y,z,vx,vy,vz]'
%%%                         [km,km,km,km/s,km/s,km/s]'
%%%


%% Parse State of Spacecraft
% Get position of spacecraft relative to the planet
Rsc_p = Xsc_p(:,1:3);

% Get velocity of spacecraft relative to the planet
Vsc_p = Xsc_p(:,4:6);


%% Get Planet State
% Set parameters to use planet location function
opts0.cBody = "Sun";
opts0.tBody = planet;

% Convert to Julian Days
tJD = SEC2JULIAN(t);

% Find planet state
Xp_s = PLANETLOC(tJD, opts0)';

%% Parse planet state
% Get position of planet relative to the sun
Rp_s = Xp_s(:,1:3);

% Get velocity of planet relative to the sun
Vp_s = Xp_s(:,4:6);

%% Calculate the State of the Spacecraft Relative to the Sun
% Calculate the position of the spacecraft relative to the sun
Rsc_s = Rsc_p + Rp_s;

% Calculate the velocity of the spacecraft relative to the sun
Vsc_s = Vsc_p + Vp_s;

% Build state vector
Xsc_s = [Rsc_s,Vsc_s];




end