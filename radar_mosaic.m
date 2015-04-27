clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% radar_mosaic.m
%
% Bonnie R. Brown, University of Hawaii Manoa, Jan 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rpth = '/Users/brbrown/data/arthur2014/radar/';
r1 = 'KLTX';    r2 = 'KMHX';
dt = '20140703';
levs = 56; x = 399;

names1 = get_radar_filenames(r1,dt,rpth,'grid');
names2 = get_radar_filenames(r2,dt,rpth,'grid');

% manual adjustment to get volumes within about 2 min of eachother
names1 = names1([2:21 23:end],:);
names2 = names2(1:end-1,:);


ref_mosaic = NaN*ones(x,x,levs,length(names1));
zdr_mosaic = ref_mosaic;
test = NaN*ones(x,x,2);
zdrtmp = NaN*ones(x,x);
ztest = NaN*ones(x,x,2);
length(names1)

for i = 1:length(names1)
    disp(['begin reading files for i = ' num2str(i)])
    ref1 = ncread([rpth r1 '/output/',dt,'/new_grid/',dt,'/',names1(i,:)],'REF');
    ref2 = ncread([rpth r2 '/output/',dt,'/new_grid/',dt,'/',names2(i,:)],'REF');
    zdr1 = ncread([rpth r1 '/output/',dt,'/new_grid/',dt,'/',names1(i,:)],'ZDR');
    zdr2 = ncread([rpth r2 '/output/',dt,'/new_grid/',dt,'/',names2(i,:)],'ZDR');
    disp(['files read, entering level loop'])
    for l = 1:56
        test(:,:,1) = squeeze(ref1(:,:,l));
        test(:,:,2) = squeeze(ref2(:,:,l));
        ztest(:,:,1) = squeeze(zdr1(:,:,l));
        ztest(:,:,2) = squeeze(zdr2(:,:,l));
        disp(['test file created'])
        
        %choose point with max reflectivity 
        [mos,ind] = max(test,[],3);
        ref_mosaic(:,:,l,i) = mos;
        disp(['max found'])
        
        for ii = 1:x
            for jj = 1:x
                if ind(ii,jj) == 1
                    zdr_mosaic(ii,jj,l,i) = ztest(ii,jj,1);
                elseif ind(ii,jj) == 2
                    zdr_mosaic(ii,jj,l,i) = ztest(ii,jj,2);
                else
                    zdr_mosaic(ii,jj,l,i) = NaN;
                end
            end
        end
        disp(['zdr retrieved'])
        
        figure(1), clf
        subplot(1,2,1)
        pcolor(mos'), shading flat
        title(['i = ' num2str(i) ' and lev = ' num2str(l)])
        subplot(1,2,2)
        pcolor(squeeze(zdr_mosaic(:,:,l,i))'), shading flat
        pause
        
        % reset tmp variabls
        test = NaN*ones(x,x,2);
        ztest = NaN*ones(x,x,2);
    end
end

newfile= 'radar_mosaic_3.nc';
% get lat and lon to write new netcdf file
lat0 = ncread([rpth r1 '/output/',dt,'/',dt,'/',names1(i,:)],'lat0');
lon0 = ncread([rpth r1 '/output/',dt,'/',dt,'/',names1(i,:)],'lon0');

nccreate(newfile,'REF','Dimensions',{'y',x,'x',x,'levels',56,'time',inf});
nccreate(newfile,'ZDR','Dimensions',{'y',x,'x',x,'levels',56,'time',inf});
nccreate(newfile,'lat0','Dimensions',{'y',x,'x',x});
nccreate(newfile,'lon0','Dimensions',{'y',x,'x',x});
ncwrite(newfile,'REF',ref_mosaic);
ncwrite(newfile,'ZDR',zdr_mosaic);
ncwrite(newfile,'lat0',lat0);
ncwrite(newfile,'lon0',lon0);


