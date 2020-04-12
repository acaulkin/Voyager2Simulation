function soiR = SOI(planet)
%%% This function returns the SOI radius in km for the planet given.
%%%
%%% INPUTS: NAME        DESCRIPTION                     SIZE        TYPE
%%%         planet  = String for desired planet SOI     NA          string
%%%                   ("Earth","Jupiter", etc.)
%%%
%%%
%%%
%%% OUTPUTS: NAME        DESCRIPTION                    SIZE        TYPE
%%%          soiR   = Radius of SOI [km]                1x1         double    
%%%
%%%

%% Calculate and return corresponding planet's SOI
switch planet
    case "Earth"
        soiR = 925000;
    case "Jupiter"
        soiR = 48200000;
    case "Saturn"
        soiR = 54800000;
    case "Uranus"
        soiR = 51800000;
    case "Neptune"
        soiR = 86600000;
end
end