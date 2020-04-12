function X = PLANETLOC(t0, opts0)
%%% This function is used to find the location of a planet(target body)
%%% relative to another body (the central body).
%%% 
%%% INPUTS: NAME        DESCRIPTION                       SIZE        TYPE
%%%         t0      = Initial time in julian days         1x1         double
%%%         opts0   = Options struct with fields
%%%                   defined below.
%%%         .cBody  = Name of central body or relative to NA          string
%%%                   ("Sun","Earth","Jupiter", etc.)
%%%         .tBody  = Name of target body                 NA          string
%%%                   ("Sun","Earth","Jupiter", etc.)
%%%
%%% OUTPUTS: NAME        DESCRIPTION                      SIZE        TYPE
%%%          X      = State vector of target planet       6x1         double
%%%                   relative to the central body                       
%%%                   [km,km,km,km/s,km/s,km/s]
%%%
%%%

%% Check Input Variables 
if isfield(opts0,'cBody')
    center = opts0.cBody;
end
if isfield(opts0,'tBody')
    target = opts0.tBody;
end

%% Calculate Desired Planet State
[pos, vel] = planetEphemeris(t0, center, target);
X = [pos, vel]';

end