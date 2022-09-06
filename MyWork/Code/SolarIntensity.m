% Function to calculate value of solar intensity for position relative to
% the Sun. Radius from Sun input taken in m.
% Code by Tanmay | 220425

function I = SolarIntensity(r)

solar_luminosity = 3.828*10^26; % https://en.wikipedia.org/wiki/Solar_luminosity

I = solar_luminosity./(4*pi*r.^2);
end