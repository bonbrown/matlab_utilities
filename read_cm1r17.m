clear all; 
%close all;
%clc;

ppth = '/Users/brbrown/data/cm1r17/nightonly/cm1out.nc';
fnum = 1;

%
% Notes: time step is one hour/60 min
%
km = 1000; % conversion factor
c0 = [0.15 0.15 0.15];

% note to self : dimensions go X x Y x Z x time (Y dimension is just 1 because axisymmetric, X is radius)
%           current: (192,1,59,166)

height = ncread(ppth,'z'); % units: km; scalar pts (59,1)
radius = ncread(ppth,'xh');  % units: km; scalar pts (192,1)
theta = ncread(ppth,'th');  % potential temperature (K)
qv = ncread(ppth,'qv');     % water vapor mixing ratio (kg/kg)
pres = ncread(ppth,'prs');  % pressure (Pa)
vt = ncread(ppth,'vinterp'); % tangential wind, cyclonic positive, interp. to scalar pts (m/s)
vr = ncread(ppth,'uinterp'); % radial wind
thetae = equiv_theta(theta,pres,qv);    % equivalent potential temperature (K) Stull (1988) approximation

[x,y,z,t] = size(theta);
to = 3600*[1:t];

% RMW at 2km height (level 12 for nightonly)
[maxv,rmi] = max(squeeze(vt(:,:,12,:)));
rmw = radius(rmi);


% cold-point tropopause
% spatial average profile
cp = 1004;      % specific heat J/(K kg)
Lv = 2260e3;    % latent heat of vaporization J/kg
k = 0.286;      % R'/cp
poo = 100000;   % reference pressure 1000 hPa 

T = theta.*(pres./poo).^k;
for ti = 1:t
    prom = squeeze(T(rmi(ti),:,:,ti));
    promT(:,ti) = prom;
    promth(:,ti) = squeeze(thetae(rmi(ti),:,:,ti));
    [tropt(ti),tropi] = min(prom);
    troph(ti) = height(tropi);
    tropth(ti) = squeeze(theta(rmi(ti),:,tropi,ti));
    tropthe(ti) = squeeze(thetae(rmi(ti),:,tropi,ti));
   
end

figure(2), hold on
plot(to/3600,troph,'r','linewidth',3);
get(gca)
set(gca,'linewidth',2,'box','on','XLim',[96 166]);
xlabel('time (hr)');
ylabel('theta (K)');
title('Cold-point tropopause temperature at RMW_{2km}')


%{
% some time averaging
te2 = squeeze(mean(thetae(:,:,:,end-60:end),4));
te1 = squeeze(mean(thetae(:,:,:,end-120:end-60),4));

vt2 = squeeze(mean(vt(:,:,:,end-60:end),4));
vt1 = squeeze(mean(vt(:,:,:,end-120:end-60),4));

vr2 = squeeze(mean(vr(:,:,:,end-60:end),4));
vr1 = squeeze(mean(vr(:,:,:,end-120:end-60),4));
%}
mp = mean(promth(:,96:166),2);
mprt = squeeze(thetae(:,:,35,96:166));
mpr = mean(mprt,2);
