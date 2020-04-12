function mu = GRAVPARAM(body)

switch body
    case "Sun"
        mu = 132712000000;
    case "Earth"
        mu = 398600;
    case "Jupiter"
        mu = 126686000;
    case "Saturn"
        mu = 37931000;
    case "Uranus"
        mu = 5794000;
    case "Neptune"
        mu = 6835100;
end

end