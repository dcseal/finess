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


fids = fopen('exact_soln.dat', 'r' );
if(fids == -1)
  disp(['File  ', 'exact_soln.dat','  not found.']);
else
    time = fscanf(fids, '%e',1);
    qtmp = fscanf(fids,'%e',[1,inf]);
    mx_ex = (length( qtmp ))/5;
    fclose(fids);
    qtmp = transpose(qtmp);
    qex = reshape(qtmp, mx_ex, 5, 1);
    dxx = 2/mx_ex;
    xex = -1+dxx/2:dxx:1';

    pex = (gamma-1)*(qex(:,5)-0.5*(qex(:,2).^2)./qex(:,1));
end    


% Density
figure(1);
clf;
pz=plot(xc,qsoln(:,1),'bo');
set(pz,'linewidth',1);
set(pz,'markersize',6);
hold off;
axis on; box on; grid off;
axis([-1 1 0 7.5]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-1:0.25:1);
set(gca,'ytick',0:1:7);
set(gca,'fontsize',16);
t1 = title(['Density at t = ',num2str(time),'     [FINESS]']); 
set(t1,'fontsize',16);
if(fids ~= -1)
    hold on;
    plot( xex, qex(:,1), '-r' );
    hold off;
end


% Pressure
figure(2);
clf;
press = (gamma-1).*(qsoln(:,5)-0.5*(qsoln(:,2).^2 + ...
                                    qsoln(:,3).^2 + qsoln(:, ...
                                                  4).^2)./qsoln(:,1));
pz=plot(xc,press,'bo');
set(pz,'linewidth',1);
set(pz,'markersize',6);
hold off;
axis on; box on; grid off;
axis([-1 1 0 0.25]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',0:0.05:0.4);
set(gca,'fontsize',16);
t1 = title(['Pressure at t = ',num2str(time),'     [FINESS]']); 
set(t1,'fontsize',16);
if(fids ~= -1)
    hold on;
    plot( xex, pex, '-r' );
    hold off;
end


% Velocity (u1)
figure(3);
clf;
pz=plot(xc,qsoln(:,2)./qsoln(:,1),'bo');
set(pz,'markersize',6);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
axis([-1 1 -1.1 1.1]);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.2:2);
set(gca,'fontsize',16);
t1 = title(['u^1(x,t) at t = ',num2str(time),'     [FINESS]']); 
set(t1,'fontsize',16);
if(fids ~= -1)
    hold on;
    plot( xex, qex(:,2)./qex(:,1), '-r' );
    hold off;
end

% Print some pretty pictures!
%figure(1);
%name1 = strcat( 'EulerTDRK4_Density', num2str(time), '.eps' );
%print(1, '-depsc2', name1 );
