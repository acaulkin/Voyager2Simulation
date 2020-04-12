%Andrew Caulkins
%AERO 343
%Voyager Mission Simulation
%12/05/2019
%__________________________________________________________________________


%% Project Statement/Summary:
%The purpose of this project is to completely simulate the voyager mission
%from 1970. It will do this by using 9 helper functions I have created that
%help define different coordinae systems/changes, sphere of influences,
%gravitational parameters, julian date, planet ephermeris data, etc...
% The main objectives of this program is to:
%
% 1. Simulate the entire Voyager 2 mission trajectory from near Earth (Aug.
% 23, 1977 7:29:11.00 UTC) to post Neptune Fly-By (Jan. 1, 1990 00:00:00 UTC).

% 2. Create four sets of plots:
% a. A full solar system perspective of the trajectory of Voyager 2 in the ICRF
% Equatorial Frame. This should include the pertinent planets orbits over the
% duration of the mission.
% b. The same as above but in a frame that is centered at Earth. This will be similar to
% one of your past homework questions, so you should be able to reuse your code.
% c. A plot showing velocity magnitude of Voyager 2 as a function of distance from
% the sun. This plot should also show the escape velocity from the Sun as a function
% of distance.
% d. Plots of each gravity assist centered around the correct central body.
% 
%
% To complete these objectives, this script will be built in 9 different sections: 
% 1. Earth to Jupiter Segment
% 2. Jupiter Fly-By
% 3. Jupiter to Saturn Segment
% 4. Saturn Fly-By
% 5. Saturn to Uranus Segment
% 6. Uranus Fly-By
% 7. Uranus to Neptune Segment
% 8. Neptune Fly-By
% 9. Post Neptune Segment
% 
% 
% This program will need to calculate the delta_v value gained due to each
% planetary flyby, and also determine how much delta_v needed to use for
% modeling/simulation errors. Using these, this code will produce an
% accurate estimation of exactly when and where Voyager 2 reached escape
% velocity from the solar system.
% 
% 
% 




clc 
clear 
close all



%% Declaration of Variables and Important Data:

%Setting options structures for the various planets:
opts0_Earth_Planet.cBody = 'Sun';
opts0_Earth_Planet.tBody = 'Earth';

opts0_Jupiter_Planet.cBody = 'Sun';
opts0_Jupiter_Planet.tBody = 'Jupiter';

opts0_Saturn_Planet.cBody = 'Sun';
opts0_Saturn_Planet.tBody = 'Saturn';

opts0_Uranus_Planet.cBody = 'Sun';
opts0_Uranus_Planet.tBody = 'Uranus';

opts0_Neptune_Planet.cBody = 'Sun';
opts0_Neptune_Planet.tBody = 'Neptune';


opts0_Jupiter.cBody = 'Jupiter';
opts0_Jupiter.tBody = 'Saturn';
opts0_Jupiter.stopCond = 1;

opts0_Jupiter_SUN.cBody = 'Sun';
opts0_Jupiter_SUN.tBody = 'Jupiter';
opts0_Jupiter_SUN.stopCond = 1;


opts0_Saturn.cBody = 'Sun';
opts0_Saturn.tBody = 'Saturn';
opts0_Saturn.stopCond = 1;

opts0_Saturn_center.cBody = 'Saturn';
opts0_Saturn_center.tBody = 'Uranus';
opts0_Saturn_center.stopCond = 1;


opts0_Uranus.cBody = 'Sun';
opts0_Uranus.tBody = 'Uranus';
opts0_Uranus.stopCond = 1;

opts0_Uranus_center.cBody = 'Uranus';
opts0_Uranus_center.tBody = 'Neptune';
opts0_Uranus_center.stopCond = 1;


opts0_Neptune.cBody = 'Sun';
opts0_Neptune.tBody = 'Neptune';
opts0_Neptune.stopCond = 1;


opts0_Neptune_center.cBody = 'Neptune';
opts0_Neptune_center.stopCond = 1;


%Importing Reference States and Seperating data into arrays:
Ref = readtable('reference_states.csv');
Ref = table2array(Ref);
Dates = Ref(:,1);
X_pos = Ref(:,2);
Y_pos = Ref(:,3);
Z_pos = Ref(:,4);
POS = [X_pos,Y_pos,Z_pos];
X_vel = Ref(:,5);
Y_vel = Ref(:,6);
Z_vel = Ref(:,7);
VEL = [X_vel,Y_vel,Z_vel];


%% Earth to Jupiter:
%Calculaing time of maneuver:
DT_1 = Dates(2) - Dates(1);
DT_1 = DT_1*86400; %Seconds

%Determining Burn Needed for re-correction:
[V1_Needed_Jupiter,V2_Needed_Jupiter] = LAMBERTSOLVE(POS(1,:),POS(2,:),DT_1,opts0_Jupiter_SUN);

%Propogating Spacecraft:
[T_1,STATE_1] = SCPROP(0,1e10,[POS(1,:),V1_Needed_Jupiter],opts0_Jupiter_SUN);


%% Jupiter Fly-by:

%Transforming from Heliocentric to Jupitercentric coordinate system:
Jupiter_State_1 = SUN2PLANET(T_1(end),STATE_1(end,:),'Jupiter');
Jupiter_State_2 = SUN2PLANET(T_1(end),STATE_1(end,:),'Jupiter');
%Propogating:
[T_2,STATE_2] = SCPROP(T_1(end),T_1(end) + 1200*86400,Jupiter_State_1',opts0_Jupiter);
% Converting: 
STATE_NEW = PLANET2SUN(T_2,STATE_2,'Jupiter');



% %The delta_V imparted by the planet is just the difference in the
% %hyperbolic excess velocities relative to Jupiter:
% VEL_START_Jupiter = State_Jupiter(4:6);
% VEL_END_Jupiter = STATE_2(end,4:6);
% 
% Delta_V_vec_Jupiter = VEL_END_Jupiter' - VEL_START_Jupiter;
% Delta_V_Jupiter = norm(Delta_V_vec_Jupiter);
% 


%% Jupiter to Saturn:

%Converting end state vector of the jupiter assist to a heliocentric state
%vector:
Saturn_State_1 = PLANET2SUN(T_2(end),STATE_2(end,:),'Jupiter');
Saturn_State_2 = PLANET2SUN(T_2(end),STATE_2(end,:),'Jupiter');
%Propogating to Saturn:
%Need to apply a correctional burn in order to head in the direction of
%Saturn. The way one could do this is through the use of Lambert's
%Algorithm because we know the heliocentric radius vector at closest 
%approach to Saturn, as well as the heliocentric radius vector at the
%beginning of our wanted trajectory. Once we know V1, we can add it to our
%current state vector's velocity and then propogate our spacecraft on its
%trajectory towards Saturn:

%Determining time between Jupiter fly-by and closest point of Saturn:
T_needed_Saturn = (Dates(3) - SEC2JULIAN(T_2(end)))*86400;
[V1_Needed_Saturn,V2_Needed_Saturn] = LAMBERTSOLVE(Saturn_State_1(1:3),POS(3,:),T_needed_Saturn,opts0_Saturn);

%Applying Correction Burn:
Saturn_State_1(4:6) = V1_Needed_Saturn;

%Proporgating to Saturn:
[T_3,STATE_3] = SCPROP(T_2(end),1e10,Saturn_State_1,opts0_Saturn);


%% Saturn Fly-by:

%Transforming from Heliocentric to Saturncentric coordinate system:
State_Saturn = SUN2PLANET(T_3(end),STATE_3(end,:),'Saturn');
%Propogating:
[T_4,STATE_4] = SCPROP(T_3(end),1e10,State_Saturn,opts0_Saturn_center);
% Converting: 
STATE_NEW_2 = PLANET2SUN(T_4,STATE_4,'Saturn');


%% Saturn to Uranus:
Uranus_State_1 = PLANET2SUN(T_4(end),STATE_4(end,:),'Saturn');
Uranus_State_2 = PLANET2SUN(T_4(end),STATE_4(end,:),'Saturn');
%Propogating to Saturn:
%Need to apply a correctional burn in order to head in the direction of
%Saturn. The way one could do this is through the use of Lambert's
%Algorithm because we know the heliocentric radius vector at closest 
%approach to Saturn, as well as the heliocentric radius vector at the
%beginning of our wanted trajectory. Once we know V1, we can add it to our
%current state vector's velocity and then propogate our spacecraft on its
%trajectory towards Saturn:

%Determining time between Saturn fly-by and closest point of Uranus:
T_needed_Uranus  = (Dates(4) - SEC2JULIAN(T_4(end)))*86400;
[V1_Needed_Uranus,V2_Needed_Uranus] = LAMBERTSOLVE(Uranus_State_1(1:3),POS(4,:),T_needed_Uranus,opts0_Uranus);

%Applying Correction Burn:
Uranus_State_1(4:6) = V1_Needed_Uranus;

%Propogating to Uranus:
[T_5,STATE_5] = SCPROP(T_4(end),1e10,Uranus_State_1,opts0_Uranus);

%% Uranus Fly-by:
%Transforming from Heliocentric to Uranuscentric coordinate system:
State_Uranus = SUN2PLANET(T_5(end),STATE_5(end,:),'Uranus');
%Propogating:
[T_6,STATE_6] = SCPROP(T_5(end),1e10,State_Uranus,opts0_Uranus_center);
% Converting: 
STATE_NEW_3 = PLANET2SUN(T_6,STATE_6,'Uranus');

%% Uranus to Neptune:
Neptune_State_1 = PLANET2SUN(T_6(end),STATE_6(end,:),'Uranus');
Neptune_State_2 = PLANET2SUN(T_6(end),STATE_6(end,:),'Uranus');
%Propogating to Saturn:
%Need to apply a correctional burn in order to head in the direction of
%Saturn. The way one could do this is through the use of Lambert's
%Algorithm because we know the heliocentric radius vector at closest 
%approach to Saturn, as well as the heliocentric radius vector at the
%beginning of our wanted trajectory. Once we know V1, we can add it to our
%current state vector's velocity and then propogate our spacecraft on its
%trajectory towards Saturn:

%Determining time between Saturn fly-by and closest point of Uranus:
T_needed_Neptune  = (Dates(5) - SEC2JULIAN(T_6(end)))*86400;
[V1_Needed_Neptune,V2_Needed_Neptune] = LAMBERTSOLVE(Neptune_State_1(1:3),POS(5,:),T_needed_Neptune,opts0_Neptune);

%Applying Correction Burn:
Neptune_State_1(4:6) = V1_Needed_Neptune;

%Propogating to Neptune:
[T_7,STATE_7] = SCPROP(T_6(end),1e10,Neptune_State_1,opts0_Neptune);

%% Neptune Fly-by:
%Transforming from Heliocentric to Neptunecentric coordinate system:
State_Neptune = SUN2PLANET(T_7(end),STATE_7(end,:),'Neptune');
%Propogating:
[T_8,STATE_8] = SCPROP(T_7(end),1e10,State_Neptune,opts0_Neptune_center);
% Converting: 
STATE_NEW_4 = PLANET2SUN(T_8,STATE_8,'Neptune');



%% Post Neptune and escape from Solar System:


%% Computing Velocity Magnitude as a Function of Distance from the Sun:
%Creating array with all state values of spacecraft in heliocentric frame:
STATES = [STATE_1;STATE_NEW;STATE_3;STATE_NEW_2;STATE_5;STATE_NEW_3;
    STATE_7;STATE_NEW_4];

%Extracting Velocity Values:
VELS = STATES(:,4:6);

%Computing the norm of the velocty values:
for i=1:size(STATES,1)
    SPEED(i) = norm(STATES(i,4:6));
end

%Computing Distance to sun:
for i=1:size(STATES,1)
    DIST(i) = norm(STATES(i,1:3));
end


%% Computing Delta_V's gained by each gravity assist and Correction Burns:
%The delta_v gained by each assist is just the difference in the hyperbolic
%excess velocities of the incoming and outbound portions of the assist
%hyperbola:

%Converting to Heliocentric


% % 
% DELTAV_Jupiter = abs(norm(STATE_2(end,4:6)) - norm(STATE_2(1,4:6)))
% DELTAV_Saturn = abs(norm(STATE_4(end,4:6)) - norm(STATE_4(1,4:6)))
% DELTAV_Uranus = abs(norm(STATE_6(end,4:6)) - norm(STATE_6(1,4:6)))
% DELTAV_Neptune = abs(norm(STATE_8(end,4:6)) - norm(STATE_8(1,4:6)))

% 


%Delta v's needed for correction burns:
DELTAV_Jupiter_Correction = abs(norm(V1_Needed_Jupiter) - norm(VEL(1,:)))
DELTAV_Saturn_Correction = abs(norm(V1_Needed_Saturn) - norm(Saturn_State_2(4:6)))
DELTAV_Uranus_Correction = abs(norm(V1_Needed_Uranus) - norm(Uranus_State_2(4:6)))
DELTAV_Neptune_Correction = abs(norm(V1_Needed_Neptune) - norm(Neptune_State_2(4:6)))



%Delta V's from gravity assist:
DELTAV_Jupiter = abs(norm(STATE_NEW(end,4:6)) - norm(STATE_NEW(1,4:6)))
DELTAV_Saturn = abs(norm(STATE_NEW_2(end,4:6)) - norm(STATE_NEW_2(1,4:6)))
DELTAV_Uranus = abs(norm(STATE_NEW_3(end,4:6)) - norm(STATE_NEW_3(1,4:6)))
DELTAV_Neptune = abs(norm(STATE_NEW_4(end,4:6)) - norm(STATE_NEW_4(1,4:6)))
DELTAV_Total = DELTAV_Jupiter + DELTAV_Saturn + DELTAV_Uranus + DELTAV_Neptune



%Calculating Escape Velocity:
V_esc = sqrt(2*GRAVPARAM('Sun')./DIST);

%% Converting to Earth Centered Frame: (Comment out if code is running too slow)

%Computing planetary times:
T_Earth = SEC2JULIAN(T_1);
T_Jupiter = SEC2JULIAN([T_1;T_2;T_3;T_4;T_5;T_6;T_7;T_8]);
T_Saturn = SEC2JULIAN([T_1;T_2;T_3;T_4;T_5;T_6;T_7;T_8]);
T_Uranus = SEC2JULIAN([T_1;T_2;T_3;T_4;T_5;T_6;T_7;T_8]);
T_Neptune = SEC2JULIAN([T_1;T_2;T_3;T_4;T_5;T_6;T_7;T_8]);

%Getting State vectors of the different planets:
X_Earth = PLANETLOC(T_Earth,opts0_Earth_Planet)';
X_Jupiter = PLANETLOC(T_Jupiter,opts0_Jupiter_Planet)';
X_Saturn = PLANETLOC(T_Saturn,opts0_Saturn_Planet)';
X_Uranus = PLANETLOC(T_Uranus,opts0_Uranus_Planet)';
X_Neptune = PLANETLOC(T_Neptune,opts0_Neptune_Planet)';

%Computing Earth reference frame based upon heliocentric planetary
%locations:
Earth_States = SUN2PLANET([T_1;T_2;T_3;T_4;T_5;T_6;T_7;T_8],STATES,'Earth');
Earth_State_Jupiter = SUN2PLANET([T_1;T_2;T_3;T_4;T_5;T_6;T_7;T_8],X_Jupiter,'Earth');
Earth_State_Saturn = SUN2PLANET([T_1;T_2;T_3;T_4;T_5;T_6;T_7;T_8],X_Saturn,'Earth');
Earth_State_Uranus = SUN2PLANET([T_1;T_2;T_3;T_4;T_5;T_6;T_7;T_8],X_Uranus,'Earth');
Earth_State_Neptune = SUN2PLANET([T_1;T_2;T_3;T_4;T_5;T_6;T_7;T_8],X_Neptune,'Earth');

%% Plotting Our Analysis:
figure(1)
%Plotting Spacecraft Trajectory:
p = plot3(STATES(:,1),STATES(:,2),STATES(:,3),'-k');
p.LineWidth = 2.0;
hold on

%Plotting Earth:
plot3(X_Earth(:,1),X_Earth(:,2),X_Earth(:,3))
%Plotting Jupiter:
plot3(X_Jupiter(:,1),X_Jupiter(:,2),X_Jupiter(:,3))
%Plotting Saturn:
plot3(X_Saturn(:,1),X_Saturn(:,2),X_Saturn(:,3))
%Plotting Uranus:
plot3(X_Uranus(:,1),X_Uranus(:,2),X_Uranus(:,3))
%Plotting Neptune:
plot3(X_Neptune(:,1),X_Neptune(:,2),X_Neptune(:,3))
%Sun:
plot3(0,0,0,'r*')

axis equal
legend('Spacecraft','Earth','Jupiter','Saturn','Uranus','Neptune','Sun')
title('Voyager 2 Mission Trajectory')


%Plotting the hyperbolic assist trajectories:
figure(2)
subplot(2,2,1)
J = plot3(STATE_2(:,1),STATE_2(:,2),STATE_2(:,3),'-k');
J.LineWidth = 2.0;
hold on
plot3(0,0,0,'r*')
title('Jupiter Fly-By')
legend('Spacecraft','Jupiter')

subplot(2,2,2)
S = plot3(STATE_4(:,1),STATE_4(:,2),STATE_4(:,3),'-k');
S.LineWidth = 2.0;
hold on
plot3(0,0,0,'r*')
title('Saturn Fly-By')
legend('Spacecraft','Saturn')

subplot(2,2,3)
U = plot3(STATE_6(:,1),STATE_6(:,2),STATE_6(:,3),'-k');
U.LineWidth = 2.0;
hold on
plot3(0,0,0,'r*')
title('Uranus Fly-By')
legend('Spacecraft','Uranus')

subplot(2,2,4)
N = plot3(STATE_8(:,1),STATE_8(:,2),STATE_8(:,3),'-k');
N.LineWidth = 2.0;
hold on
plot3(0,0,0,'r*')
title('Neptune Fly-By')
legend('Spacecraft','Neptune')

%Plotting the magnitude of velocity as a function of distance from the sun
%as well as the escape velocity of the solar system based on distance:
figure(3)
L = plot(DIST,SPEED,'-k');
hold on
plot(DIST,V_esc,'--r')
title('Velocity as a Function of Distance from the Sun')
xlabel('Distance from Sun [km]')
ylabel('Velocity Magnitude [km/s]')
legend('Spacecraft Velocity','Solar System Escape Velocity')

%Plotting Earth centered frame analysis:
figure(4)
O = plot3(Earth_States(:,1),Earth_States(:,2),Earth_States(:,3),'-k');
O.LineWidth = 2.0;
hold on
plot3(Earth_State_Jupiter(:,1),Earth_State_Jupiter(:,2),Earth_State_Jupiter(:,3));
plot3(Earth_State_Saturn(:,1),Earth_State_Saturn(:,2),Earth_State_Saturn(:,3));
plot3(Earth_State_Uranus(:,1),Earth_State_Uranus(:,2),Earth_State_Uranus(:,3));
plot3(Earth_State_Neptune(:,1),Earth_State_Neptune(:,2),Earth_State_Neptune(:,3));
plot3(0,0,0,'r*')

axis equal
legend('Spacecraft','Jupiter','Saturn','Uranus','Neptune','Earth')
title('Voyager 2 Mission Trajectory')




