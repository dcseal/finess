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
ion_mass   = sscanf(INI.plasma.ion_mass, '%e');
elc_mass   = ion_mass / mass_ratio;

cs_light   = sscanf(INI.maxwell.cs_light, '%e');

fprintf(1,'    mass ratio = %f  ; speed of light = %f \n', mass_ratio, cs_light );


figure(1);
clf;
pz=plot(xc,qsoln(:,1),'bo');
hold on;
plot(xc,qsoln(:,6),'rx');
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
