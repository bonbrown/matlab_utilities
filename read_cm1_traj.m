clear all
%close all
clc
%
% Notes: time step is 60 sec
%   parcels are released at 4 days (345600 sec) and go until ~6.9 days
%   (596520 sec)
%
ppth = '/Users/brbrown/data/cm1r17/nightonly/cm1out_pdata.nc';
fnum = 2;

km = 1000; % conversion factor
%c0 = [0 0.15 1.0];
%c00 = [1.0 0.15 0];
c0 = [0.15 0.15 0.15];

height = ncread(ppth,'z'); % units: m; dimensions: parcel # x time (600 x 4183)
radius = ncread(ppth,'x');  % units: m
theta = ncread(ppth,'th');  % potential temperature (K)
qv = ncread(ppth,'qv');     % water vapor mixing ratio (kg/kg)
pres = ncread(ppth,'prs');  % pressure (Pa)
u = ncread(ppth,'u');   % radial wind (m/s)
v = ncread(ppth,'v');   % tangential wind (m/s)
thetae = equiv_theta(theta,pres,qv);    % equivalent potential temperature (K) Stull (1988) approximation

[p,t] = size(theta);
tp = 60*[1:t];

c = 0;

for i = 1:1:600 % [1:30 100:130 200:230 300:330 400:430 500:530 600]
    
    c = c + 3;
    
    te = squeeze(thetae(i,:));
    vt = squeeze(v(i,:));
    ur = squeeze(u(i,:));
    r = squeeze(radius(i,:))/km;
    h = squeeze(height(i,:))/km;
    
   %mte = max(thetae(i,:));
    %mv = min(u(i,:));
    
    i0 = find(vt<=0,1,'first');
    tet = h(i0);
    
    figure(fnum), hold on
    
    %subplot(2,2,4)
    %plot(vt,te,'color',[0 c0(2)*(1+c*0.025) 1],'HandleVisibility','off');
    %plot(ur,te,'color',[1 c00(2)*(1+c*0.025) 0],'HandleVisibility','off');
    %plot(h,te,'color',c0*(1+c*0.025),'HandleVisibility','off');
    
    if i==1
        plot((345600+tp(i0))/3600,tet,'r+','linewidth',3);
    else
        plot((345600+tp(i0))/3600,tet,'r+','linewidth',3,'HandleVisibility','off');
    end
    
    
end

%{
figure(fnum), hold on

xlabel('m/s');
ylabel('equiv. pot. temperature (K)');
get(gca);  
set(gca,'linewidth',2,'box','on','fontsize',10,'fontweight','bold')
%title('(d) Night-only run');
%}