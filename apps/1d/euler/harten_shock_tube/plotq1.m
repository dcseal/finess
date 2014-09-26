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
%% ----
%% HACK
%% ----
gamma = 1.4;
%% ----

% pressure
press = (gamma-1)*(qsoln(:,5)-0.5*(qsoln(:,2).^2+...
                                   qsoln(:,3).^2+...
                                   qsoln(:,4).^2)./qsoln(:,1));

% -------------------------------------------- %
% Exact solution (at the final time):
% Uncomment if not desired.  The data for this 
% file was produced by HYPERPYWS.
fids = fopen('exact_soln.dat', 'r' );
if(fids == -1)
  disp(['File  ', 'rk4-soln-100.dat','  not found.']);
else
    qex = fscanf(fids,'%e',[4,inf]);
    fclose( fids );
    mx_ex = length( qex );
    xex = qex( 1, : )';  qex = qex( 2:4, : )';

    pex = (gamma-1)*(qex(:,3)-0.5*(qex(:,2).^2)./qex(:,1));
end    
% -------------------------------------------- %

% -------------------------------------------- %
% Secondary solution to compare result to
fids_extra = fopen('rk4-soln-100.dat', 'r' );
if(fids_extra == -1)
  disp(['File  ', 'exact_soln.dat','  not found.']);
else
    qextra = fscanf(fids,'%e',[4,inf]);
    fclose( fids );
    mx_extra = length( qextra );
    x_extra = qextra( 1, : )';  qextra = qextra( 2:4, : )';
    pextra  = (gamma-1)*(qextra(:,3)-0.5*(qextra(:,2).^2)./qextra(:,1));
end    
% -------------------------------------------- %


figure(1);
clf;
pz=plot(xc,qsoln(:,1),'bo');
set(pz,'linewidth',1);
set(pz,'markersize',6);
hold off;
axis on; box on; grid off;
axis([-0.01 1.01 0.29 1.4]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.2:5);
set(gca,'ytick',-3:0.2:12);
set(gca,'fontsize',16);
%t1 = title(['Density at t = ',num2str(time),'     [FINESS]']); 
t1 = title(['Density']);
set(t1,'fontsize',16);

if(fids_extra ~= -1)
    hold on;
    plot( x_extra, qextra(:,1), 'g+' );
    hold off;
end
if(fids ~= -1)
    hold on;
    plot( xex, qex(:,1), '-r' );
    hold off;
end
l1 = legend('Lax-Wendroff', 'Runge-Kutta', 'Exact');
set(l1, 'Location', 'NorthWest');

figure(2);
clf;
pz=plot(xc, press, 'bo');
set(pz,'linewidth',1);
set(pz,'markersize',6);
hold off;
axis on; box on; grid off;
axis([-0.01 1.01  0.4 3.66]);

set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.20:5);
set(gca,'ytick',-0:0.5:40);
set(gca,'fontsize',16);
%t1 = title(['Pressure at t = ',num2str(time),'     [FINESS]']);
t1 = title(['Pressure']);
set(t1,'fontsize',16);

if(fids_extra ~= -1)
    hold on;
    plot( x_extra, pextra(:,1), 'g+' );
    hold off;
end
if(fids ~= -1)
    hold on;
    plot( xex, pex, '-r' );
    hold off;
end
l2 = legend('Lax-Wendroff', 'Runge-Kutta', 'Exact');
set(l2,'Location', 'Best');

figure(3);
clf;
pz=plot(xc, qsoln(:,2)./qsoln(:,1), 'bo');
set(pz,'markersize',6);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
axis([-0.01 1.01 -0.01 1.66]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.2:5);
set(gca,'ytick',-1:0.5:5);
set(gca,'fontsize',16);
%t1 = title(['u^1(x,t) at t = ',num2str(time),'     [FINESS]']); 
t1 = title(['Velocity']);
set(t1,'fontsize',16);


if(fids_extra ~= -1)
    hold on;
    plot( x_extra, qextra(:,2)./qextra(:,1), 'g+' );
    hold off;
end
if(fids ~= -1)
    hold on;
    plot( xex, qex(:,2)./qex(:,1), '-r' );
    hold off;
end
l3 = legend('Lax-Wendroff', 'Runge-Kutta', 'Exact');
set( l3, 'Location', 'South' );

% Save the pretty pictures!
print(1, '-depsc', 'harten_density.eps'  );
print(2, '-depsc', 'harten_pressure.eps' );
print(3, '-depsc', 'harten_velocity.eps' );
