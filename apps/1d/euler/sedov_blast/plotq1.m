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

% gas constant
if( isfield( INI.euler, 'gamma' ) )
    gamma = sscanf( INI.euler.gamma, '%f' );
else
    gamma = 1.4;
end

% Density
figure(1);
clf;
pz=plot(xc,qsoln(:,1),'bo');
set(pz,'linewidth',1);
set(pz,'markersize',8);
hold off;
axis on; box on; grid off;
axis([-2 2 0 6.5]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.5:2);
set(gca,'ytick',0:1:7);
set(gca,'fontsize',16);
t1 = title(['Density at t = ',num2str(time),'     [FINESS]']); 
set(t1,'fontsize',16);

% Pressure
figure(2);
clf;
press = (gamma-1).*(qsoln(:,5)-0.5*(qsoln(:,2).^2 + ...
                                    qsoln(:,3).^2 + qsoln(:, ...
                                                  4).^2)./qsoln(:,1));
pz=plot(xc,press,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',8);
hold off;
axis on; box on; grid off;
axis([-2 2 0 8e5]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.5:2);
set(gca,'ytick',0:1e5:10e5);
set(gca,'fontsize',16);
t1 = title(['Pressure at t = ',num2str(time),'     [FINESS]']); 
set(t1,'fontsize',16);

% Velocity (u1)
figure(3);
clf;
pz=plot(xc,qsoln(:,2)./qsoln(:,1),'bo');
set(pz,'markersize',8);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
axis([-2 2 -1000 1000]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.5:2);
set(gca,'ytick',-2000:200:2000);
set(gca,'fontsize',16);
t1 = title(['u^1(x,t) at t = ',num2str(time),'     [FINESS]']); 
set(t1,'fontsize',16);

% Print some pretty pictures!
%figure(1);
%name1 = strcat( 'EulerTDRK4_Density', num2str(time), '.eps' );
%print(1, '-depsc2', name1 );
