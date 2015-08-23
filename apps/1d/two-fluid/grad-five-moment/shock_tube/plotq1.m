%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%              mx:  number of points
%%    [xlow,xhigh]:  min/max values of grid
%%            meqn:  number of equations
%%            maux:  number of aux components
%%           meth1:  spatial order of accuracy
%%
%%   Grid information:
%%              xc: grid points (cell centers), size = (mx,my)
%%
%%   Solution information:
%%           qsoln:  solution sampled on mesh, size = (mx,meqn)
%%             aux:  aux components sampled on mesh, size = (mx,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Model parameters (and derived parameters too)
gamma      = sscanf(INI.plasma.gamma, '%e');
mass_ratio = sscanf(INI.plasma.mass_ratio, '%e');
ion_mass   = 1.0;
elc_mass   = ion_mass / mass_ratio;
cs_light   = sscanf(INI.maxwell.cs_light, '%e');

fprintf(1,'    mass ratio = %f  ; speed of light = %f \n', mass_ratio, cs_light );

% Indices for components
ind_rho_i     = 1;  %  rho (ion)
ind_M1_i      = 2;  %  1-momentum (ion)
ind_M2_i      = 3;  %  2-momentum (ion)
ind_M3_i      = 4;  %  3-momentum (ion)
ind_N_i       = 5;  %  energy (ion)
ind_rho_e     = 6;  %  rho (elc)
ind_M1_e      = 7;  %  1-momentum (elc)
ind_M2_e      = 8;  %  2-momentum (elc)
ind_M3_e      = 9;  %  3-momentum (elc)
ind_N_e       = 10; %  energy (elc)
ind_B1        = 11; %  1-magnetic field
ind_B2        = 12; %  2-magnetic field
ind_B3        = 13; %  3-magnetic field
ind_E1        = 14; %  1-electric field
ind_E2        = 15; %  2-electric field
ind_E3        = 16; %  3-electric field
ind_psi       = 17; %  B-field cleaning
ind_phi       = 18; %  E-field cleaning
ind_entropy_i = 19; %  entropy tracking
ind_entropy_e = 20; %  entropy tracking

J1_i =  qsoln(:,ind_M1_i) / ion_mass;
J1_e = -qsoln(:,ind_M1_e) / elc_mass;
J1  = J1_i + J1_e;

%   u1_i = qsoln(:,5) / q(:,1);
%   u1_i = qsoln(:,2) / ion_mass;

figure(1);
clf;
pz=plot(xc,qsoln(:,ind_rho_i),'bo');
hold on;
plot(xc,qsoln(:,ind_rho_e),'rx');
hold off;
set(pz,'linewidth',1);
set(pz,'markersize',8);
hold off;
axis on; box on; grid off;
axis([-0.1 13 0 1.1]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',0:1:15);
set(gca,'ytick',0:.1:2);
set(gca,'fontsize',16);
t1 = title(['Density at t = ',num2str(time),'     [FINESS]']); 
set(t1,'fontsize',16);
legend('ion density', 'electron density');

figure(2);
clf;
pz=plot(xc,qsoln(:,ind_E1),'bo');
%hold on;
%plot(xc,qsoln(:,6),'rx');
%hold off;
%   set(pz,'linewidth',1);
%   set(pz,'markersize',8);
%   hold off;
axis on; box on; grid off;
axis auto;
%axis([-0.1 13 0 1.1]);
%   set(gca,'plotboxaspectratio',[1.5 1 1]);
%   set(gca,'xtick',0:1:15);
%   set(gca,'ytick',0:.1:2);
%   set(gca,'fontsize',16);
%   t1 = title(['Density at t = ',num2str(time),'     [FINESS]']); 
%   set(t1,'fontsize',16);
%   legend('ion density', 'electron density');

figure(3);
clf;
pz=plot(xc, J1, 'bo');
set(pz,'linewidth',1);
hold on;
set(pz,'markersize',8);
plot(xc, J1_i, 'rx' );
hold on;
plot(xc, J1_e, 'go' );
hold off;
axis on; box on; grid off;
%axis([-0.1 13 0 1.1]);
axis auto;
%set(gca,'plotboxaspectratio',[1.5 1 1]);
%set(gca,'xtick',0:1:15);
%set(gca,'ytick',0:.1:2);
set(gca,'fontsize',16);
t1 = title(['Current at t = ',num2str(time),'     [FINESS]']); 
set(t1,'fontsize',16);


%   figure(2);
%   clf;
%   press = (gamma-1).*(qsoln(:,5)-0.5*(qsoln(:,2).^2 + ...
%                                       qsoln(:,3).^2 + qsoln(:, ...
%                                                     4).^2)./qsoln(:,1));
%   pz=plot(xc,press,'bo');
%   set(pz,'linewidth',1);
%   set(pz,'markersize',8);
%   hold off;
%   axis on; box on; grid off;
%   %axis([-5 5 0 12]);
%   set(gca,'plotboxaspectratio',[1.5 1 1]);
%   set(gca,'xtick',-5:2.5:5);
%   set(gca,'ytick',0:2:12);
%   set(gca,'fontsize',16);
%   t1 = title(['Pressure at t = ',num2str(time),'     [FINESS]']); 
%   set(t1,'fontsize',16);
%   if(fids ~= -1)
%       hold on;
%       plot( xex, pex, '-r' );
%       hold off;
%   end

%   figure(3);
%   clf;
%   pz=plot(xc,qsoln(:,2)./qsoln(:,1),'bo');
%   set(pz,'markersize',8);
%   set(pz,'linewidth',1);
%   hold off;
%   axis on; box on; grid off;
%   %axis([-5 5 -0.5 3]);
%   set(gca,'plotboxaspectratio',[1.5 1 1]);
%   set(gca,'xtick',-5:2.5:5);
%   set(gca,'ytick',0:1:3);
%   set(gca,'fontsize',16);
%   t1 = title(['u^1(x,t) at t = ',num2str(time),'     [FINESS]']); 
%   %t1 = title(['Velocity']);
%   set(t1,'fontsize',16);

%   if(fids ~= -1)
%       hold on;
%       plot( xex, qex(:,2)./qex(:,1), '-r' );
%       hold off;
%   end

%   % Save the pretty pictures!
%   print(1, '-depsc', 'shock_entropy_density.eps'  );
%   print(2, '-depsc', 'shock_entropy_pressure.eps' );
%   print(3, '-depsc', 'shock_entropy_velocity.eps' );
