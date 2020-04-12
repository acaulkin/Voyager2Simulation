function julianDate = SEC2JULIAN(t)
%%% This function returns the julian date calculated from a time in seconds
%%% given from the starting julian date that is defined inside of this
%%% function. This date will be the date of launch of Voyager 2.

%% Set Starting Julian Date
J0 = juliandate(1977,8,23,7,29,11);

%% Transform Seconds to Days
sec2hrs = 1/3600;
hrs2days = 1/24;
days = t*sec2hrs*hrs2days;

%% Calculate Julian Date of Time Given
julianDate = J0+days;

end