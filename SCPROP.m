function [T,X] = SCPROP(t0, tf, X0, opts0)
%%% This function is to be used for propagating the state of a spacecraft
%%% around some central body using the conditions defined by the inputs and
%%% options. This may be from an intial time to a final time or to some
%%% physical stopping point as defined in the options struct.
%%% 
%%% INPUTS: NAME        DESCRIPTION                         SIZE     TYPE
%%%              t0 = Initial time (s)                      1x1      double
%%%              tf = Final time (s)                        1x1      double
%%%              X0 = Initial state [x,y,z,vx,vy,vz]'       6x1      double
%%%                    [km,km,km,km/s,km/s,km/s]'
%%%                    In the frame of the central body
%%%           opts0 = Options struct with fields
%%%                   defined below.
%%%          .cBody = Name of central body orbited          NA       string
%%%                   ("Sun","Earth","Jupiter", etc.)
%%%          .tBody = Name of target body if needed         NA       string
%%%       .stopCond = Code stating the stop condition       1x1      int
%%%                 = 1, stop at the SOI
%%%
%%%
%%% OUTPUTS: NAME        DESCRIPTION                     SIZE        TYPE
%%%               T = Time vector history (s)            Nx1         double
%%%               X = State vector history [s]           Nx6         double
%%%
%%%


%% Gather General Constants and Parameters
if isfield(opts0,'cBody')
    mu = GRAVPARAM(opts0.cBody);
    
    if opts0.cBody~="Sun"
        soiRcentral = SOI(opts0.cBody);
    else
        soiRcentral = false;
    end
else
    error('A central body was not specified for this propagation segment')
end

if isfield(opts0,'tBody')
    soiRtarget = SOI(opts0.tBody);
end
    
if isfield(opts0,'stopCond')
    stopCode = opts0.stopCond;
end

%% Build ODE Options
% Define a general options structure
OPTIONS = odeset('AbsTol',1E-25,'RelTol',2.23E-14);

% Define event function for propagation if there is one
if isfield(opts0,'cBody') && isfield(opts0,'stopCond')
    if stopCode==1 
        if opts0.cBody=="Sun" && isfield(opts0,'tBody')
            OPTIONS.Events = @(T,X)SOI_IN_STOP(T,X,soiRtarget,opts0.tBody);                
        elseif opts0.cBody~="Sun"
            OPTIONS.Events = @(T,X)SOI_OUT_STOP(T,X,soiRcentral);
        else
            error("Correct information was not given to be able to apply stop condition")
        end
    end   
end

% Time parameters
tspan = [t0 tf];
%tspan = t0:1000:tf;
%% Run ODE113 Solver
[T,X] = ode113(@(T,X) PROP(T,X,mu),tspan,X0,OPTIONS);

%% Dynamics Function for Two-Body Problem
    function dX = PROP(T,X,mu)
        % Parse State Vector
        r  = X(1:3);                % Position vector [km]
        v  = X(4:6);                % Velocity vector [km/s]
        
        % Calculate distance from Central Body to Spacecraft
        rn = norm(r);               % Distance [km]
                   
        % Building State Time Derivative
        dX(1:3)  = v;
        dX(4:6)  = -mu/rn^3*r;  
        dX = dX(:);              % Convert to vector column format   
    end

%% Stop Condition at the SOI for Spacecraft traveling from the outside in
    function [value,isterminal,direction] = SOI_IN_STOP(T,X,soiR,body)
        %%% Integration stops when the trajectory pierces the SOI of the body
        % Transfer Coordinates Into Planet Frame
        Xsc_p = SUN2PLANET(T,X',body);
        
        % Get Distance from Planet in Planet Frame
        Rsc_p = [Xsc_p(1),Xsc_p(2),Xsc_p(3)];
        Dsc_p = norm(Rsc_p);
        
        % Set up stopping conditions
        %value = ~((Dsc_p-soiR)<0);     % Zero state value triggers event
        value = Dsc_p-soiR;
        isterminal = 1;                % Stop integration at event
        direction  = -1;    
    end

%% Stop Condition at the SOI for Spacecraft traveling from the inside out
    function [value,isterminal,direction] = SOI_OUT_STOP(T,X,soiR)
        %%% Integration stops when the trajectory pierces the SOI of the body
        % Parse Position and Distance of SC from the Body
        Rsc_p = [X(1),X(2),X(3)];
        Dsc_p = norm(Rsc_p);
        
        % Set up stopping
        value = Dsc_p-soiR;      % Zero state value triggers event
        isterminal = 1;         % Stop integration at event
        direction  = 1;    
    end
end