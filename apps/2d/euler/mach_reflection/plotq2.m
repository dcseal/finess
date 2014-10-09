%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%                  mx, my:  number of points in each coordinate direction
%% [xlow,xhigh,ylow,yhigh]:  min/max values of grid
%%                    meqn:  number of equations
%%                    maux:  number of aux components
%%
%%   Grid information:
%%       (xc,yc): grid points (cell centers), size = (mx,my)
%%
%%   Solution information:
%%         qsoln:  solution sampled on mesh, size = (mx,my,meqn)
%%           aux:  aux components sampled on mesh, size = (mx,my,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gas constant
global INI;
gamma_gas = sscanf(INI.euler.gamma, '%e');

figure(1);
clf;
pcolor(xc,yc,qsoln(:,:,m));
shading flat;
yrbcolormap
axis on; box on; grid off;
axis('equal');
axis([-0.05 3.25 -0.05 1.05]);
set(gca,'xtick',-4:0.5:4);
set(gca,'ytick',-4:0.5:4);
set(gca,'fontsize',16);
%t1 = title(['\rho(t = ', num2str(time,'%2.2e'), ', x)     [FINESS]']); 
%t1 = title(['\rho(t, x, y)     [FINESS]']); 
t1 = title(['Density at t = ',num2str(time),'     [FINESS]']); 
set(t1,'fontsize',16);
colorbar;
caxis([1,25]);
   
   figure(2)
   clf
   plot(xc(:,1),qsoln(:,25,1),'b-');

   % Grey scale figure plot that you see in a lot of papers.  For example, see 
   % Figs 4.1 and 4.2 in 
   %
   %    "The Runge–Kutta Discontinuous Galerkin Method for Conservation
   %     Laws V", Cockburn and Chi-Wang Shu, J. Comp. Phys., 141, 199–224 (1998).
   %  
   figure(3);
   clf;
   %contour(xc, yc, qsoln(:,:,m), linspace(1.3965, 22.682,30), '-k' );
   % Woodward and Collela's contour lines:
   contour(xc, yc, qsoln(:,:,m), linspace(1.728, 20.74,30), '-k' );
   axis on; box on; grid off;
   axis('equal');
   axis([-0.05 3.25 -0.05 1.05]);
   set(gca,'xtick',-4:0.5:4);
   set(gca,'ytick',-4:0.5:4);
   set(gca,'fontsize',16);
   %t1 = title(['\rho(t,x,y) at t = ',num2str(time),'     [FINESS]']); 
   %t1 = title(['Density at t = ',num2str(time),'     [FINESS]']); 
   t1 = title(['Density']);
   set(t1,'fontsize',16);
 
figure(1)

% n1 = frame number
%fname = strcat( strcat( 'density', num2str(n1, '%02d' ) ), '.jpg' );
%print(1, '-djpeg', fname );
%fname = strcat( strcat( 'density-contour', num2str(n1, '%02d' ) ), '.eps' );
%print(3, '-deps', fname  );

