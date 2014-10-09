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
caxis([0.8,10]);
   
% Grey scale figure plot that you see in a lot of papers.
figure(3);
clf;
contour(xc, yc, qsoln(:,:,m), linspace(0.2365, 5.647,30), '-k' );
contour(xc, yc, qsoln(:,:,m), linspace(0.2568, 6.607,30), '-k' );  % Jiang and Shu lines
axis on; box on; grid off;
axis('equal');
axis([-0.05 3.05 -0.05 1.05]);
set(gca,'xtick',-4:0.5:4);
set(gca,'ytick',-4:0.5:4);
set(gca,'fontsize',16);
%t1 = title(['\rho(t,x,y) at t = ',num2str(time),'     [FINESS]']); 
t1 = title(['Density at t = ',num2str(time),'     [FINESS]']); 
%t1 = title(['Density']);
set(t1,'fontsize',16);

figure(1)

% n1 = frame number
%fname = strcat( strcat( 'density', num2str(n1, '%02d' ) ), '.jpg' );
%print(1, '-djpeg', fname );
%fname = strcat( strcat( 'density-contour', num2str(n1, '%02d' ) ), '.eps' );
%print(3, '-deps', fname  );

