
function [wavefield_gradient,gather]=TMGenerateData(srcloc,recloc,srcpulse,T,isrc,ep,mu,sig,dt,dx,dz,x,z,outstep, plotopt,dt_wf,nt_wf)
% Example for running TM_model2d.m
% (a 2-D, FDTD, reflection GPR modeling code in MATLAB)
%
% by James Irving
% July 2005

% load the earth model example from file
%load true_model;
% calculate minimum and maximum relative permittivity and permeability
% in the model (to be used in finddx.m and finddt.m) 
epmin = min(min(ep)); % 1
epmax = max(max(ep)); % 25
mumin = min(min(mu)); % 1
mumax = max(max(mu)); % 1

% create time (s) and source pulse vectors for use with finddx.m
% (set dt very small to not alias frequency components;  set again below)
% (maximum time should be set to that whole source pulse is included)
% t=0:1e-10:100e-9;
% srcpulse=blackharrispulse(100e6,t);

% use finddx.m to determine maximum possible spatial field discretization
% (in order to avoid numerical dispersion)
% [dx,wlmin,fmax] = finddx(epmax,mumax,srcpulse,t,0.02);
% disp(' ');
% disp(['Maximum frequency contained in source pulse = ',num2str(fmax/1e6),' MHz']);
% disp(['Minimum wavelength in simulation grid = ',num2str(wlmin),' m']);
% disp(['Maximum possible electric/magnetic field discretization (dx,dz) = ',num2str(dx),' m']);
% disp(['Maximum possible electrical property discretization (dx/2,dz/2) = ',num2str(dx/2),' m']);
% disp(' ');

% set dx and dz here (m) using the above results as a guide
% dx = 0.04;
% dz = 0.04;
%disp(['Using dx = ',num2str(dx),' m, dz = ',num2str(dz),' m']); % dx=0.1m dz=0.1m

% find the maximum possible time step using this dx and dz
% (in order to avoid numerical instability)
dtmax = finddt(epmin,mumin,dx,dz);   %2.0217e-10
%disp(['Maximum possible time step with this discretization = ',num2str(dtmax/1e-9),' ns']);
%disp(' ');

% set proper dt here (s) using the above results as a guide
% dt = 8e-11;
%disp(['Using dt = ',num2str(dt/1e-9),' ns']);  % dt=0.08ns
%disp(' ');

% create time vector (s) and corresponding source pulse
% (using the proper values of dt and tmax this time)
t=0:dt:T; % 2501 st                          
% srcpulse = blackharriswave(100e6,t);    
% srcpulse = blackharrispulse(100e6,t);
% interpolate electrical property grids to proper spatial discretization
% NOTE:  we MUST use dx/2 here because we're dealing with electrical property matrices
%disp('Interpolating electrical property matrices...');
%disp(' ');
x2 = min(x):dx/2:max(x);
z2 = min(z):dx/2:max(z);
ep2 = gridinterp(ep,x,z,x2,z2,'nearest');%临近域线性插值
mu2 = gridinterp(mu,x,z,x2,z2,'nearest');
sig2 = gridinterp(sig,x,z,x2,z2,'nearest');

% plot electrical property grids to ensure that interpolation was done properly
% figure; subplot(2,1,1);
% imagesc(x,z,ep'); axis image; colorbar
% xlabel('x (m)'); ylabel('z (m)');
% title('Original \epsilon_r matrix');
% subplot(2,1,2)
% imagesc(x2,z2,ep2'); axis image; colorbar
% xlabel('x (m)'); ylabel('z (m)');
% title('Interpolated \epsilon_r matrix');
% %
% figure; subplot(2,1,1);
% imagesc(x,z,mu'); axis image; colorbar
% xlabel('x (m)'); ylabel('z (m)');
% title('Original \mu_r matrix');
% subplot(2,1,2)
% imagesc(x2,z2,mu2'); axis image; colorbar
% xlabel('x (m)'); ylabel('z (m)');
% title('Interpolated \mu_r matrix');
% %
% figure; subplot(2,1,1);
% imagesc(x,z,sig'); axis image; colorbar
% xlabel('x (m)'); ylabel('z (m)');
% title('Original \sigma matrix');
% subplot(2,1,2)
% imagesc(x2,z2,sig2'); axis image; colorbar
% xlabel('x (m)'); ylabel('z (m)');
% title('Interpolated \sigma matrix');

% pad electrical property matrices for PML absorbing boundaries
npml = 10;  % number of PML boundary layers
[ep3,x3,z3] = padgrid(ep2,x2,z2,2*npml+1);%padgrid为了扩充边界
[mu3,x3,z3] = padgrid(mu2,x2,z2,2*npml+1);
[sig3,x3,z3] = padgrid(sig2,x2,z2,2*npml+1);

% clear unnecessary matrices taking up memory
clear  x2  z2 ep ep2 mu mu2 sig sig2 

% create source and receiver location matrices
% (rows are [x location (m), z location (m)])
% srcx = (0:2:18)';
% srcz = 0*ones(size(srcx));
% recx = srcx + 1;
% recz = srcz;
% srcloc = [srcx srcz];
% recloc = [recx recz];
% 
% % set some output and plotting parameters
% outstep = 4;
% plotopt = [1 50 0.002];
% 
% % pause
% disp('Press any key to begin simulation...');
% disp(' ');
% pause;
close all

% run the simulation
% ---tic;
[wavefield_gradient,gather] = TM_model2d1(ep3,mu3,sig3,x3,z3,srcloc,recloc,srcpulse,t,npml,outstep,plotopt,x,z,dt_wf,nt_wf);
%disp(' ');
%disp(['Total running time = ',num2str(toc/3600),' hours']);

%  save(['Gather00_',num2str(isrc),'.mat'],'gather','tout','srcx','srcz','recx','recz','dt','dx','dz','x','z','-v7.3')
%  save(['Wavefield00_',num2str(isrc),'.mat'],'wavefield','tout','srcx','srcz','recx','recz','dt','dx','dz','x','z','-v7.3')
%   save(['Gather00_',num2str(isrc),'.mat'],'gather')
%   save(['Wavefield00_',num2str(isrc),'.mat'],'wavefield')
%save(['Gather00_',num2str(isrc),'.mat'],'gather',tout)
%save(['Gather00_',num2str(isrc),'.mat'],'gather','tout')
%save(['Wavefield00_',num2str(isrc),'.mat'],'wavefield','tout')
end

%% extract common offset reflection GPR data from multi-offset data cube and plot the results
% for i=1:length(srcx);
%     co_data(:,i) = gather(:,i,i);
% end
% pos = (srcx+recx)/2;
% figure; subplot(2,2,[1 2]);
% imagesc(pos,tout*1e9,co_data);
% axis([0 20 0 250]);
% set(gca,'plotboxaspectratio',[2 1 1]);
% caxis([-5e-4 5e-4]);
% colormap('gray');
% xlabel('Position (m)');
% ylabel('Time (ns)');