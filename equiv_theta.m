function equiv_theta = equiv_theta(theta,pres,mxratio)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs: potential temperature (theta; K)
%         pressure (pres; Pa)
%         water vapor mixing ratio (mxratio; kg/kg)
%
%   Approximate calculation of equivalent potential temperature from
%       Stull (1988) sec.13.1 page 546 
%
%   Bonnie R. Brown, University of Hawai'i Manoa, April 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cp = 1004;      % specific heat J/(K kg)
Lv = 2260e3;    % latent heat of vaporization J/kg
k = 0.286;      % R'/cp
poo = 100000;   % reference pressure 1000 hPa 

T = theta.*(pres./poo).^k;

equiv_theta = (T + (Lv/cp)*mxratio).*(poo./pres).^k;

