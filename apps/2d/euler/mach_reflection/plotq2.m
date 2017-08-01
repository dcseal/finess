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

% gas constant (this is a good example of how to pull items from the
% parameters.ini file)
global INI;
gamma_gas = sscanf(INI.euler.gamma, '%e');

% ---------------------------------------------------------------------------- %
% Colorplot of the variable of choice
% ---------------------------------------------------------------------------- %

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

% ---------------------------------------------------------------------------- %
% Schlieren plot for density (or any other variable)
% ---------------------------------------------------------------------------- %

% Grey scale figure plot that you see in a lot of papers.
% figure('Position', [100, 100, 1049, 895]);
figure(2);
clf;

rho=qsoln(:,:,1);

%colors=linspace(1.0,3.0,10);
h=xc(2,1)-xc(1,1);
[px,py]=gradient(rho,h,h);
v=sqrt(px.^2+py.^2)';

fprintf('Smallest and largest value of the computed gradient = %f %f\n', max(max(v)), min(min(v)) );
colormap( flipud(gray(2048)).^10 )
hp = pcolor(xc,yc,v');
set(hp, 'EdgeColor', 'none');

% Fancier labels for the plots
axis on; box on; grid off;
axis('equal');
axis([-0.05 3.25 -0.05 1.05]);
set(gca,'xtick',-4:0.5:4);
set(gca,'ytick',-4:0.5:4);
set(gca,'fontsize',16);
t1 = title(['Density at t = ',num2str(time),'     [FINESS]']); 
set(t1,'fontsize',16);

% ---------------------------------------------------------------------------- %
% Grey scale topographical figure plot (this smooths out a lot of features but
% seems to be the standard thing to plot)
% ---------------------------------------------------------------------------- %

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
t1 = title(['Density at t = ',num2str(time),'     [FINESS]']); 
%t1 = title(['Density']);
set(t1,'fontsize',16);

% ---------------------------------------------------------------------------- %
% Slice of the solution (to check for oscillations)
% ---------------------------------------------------------------------------- %

figure(4)
clf
plot( xc(:,1), qsoln(:,25,1), 'b-');

set(gca,'xtick',-4:0.5:4);
set(gca,'fontsize',16);
xlabel('x');
%ylabel('\rho', 'rot',0);
ylabel('Density');
set(gca,'plotboxaspectratio',[2 1 1]);

axis([-0.05 3.25 -0.05 20.05]);

t1 = title(['Density at y = ',num2str(xc(1,25)),'     [FINESS]']); 
set(t1,'fontsize',16);

figure(1)

% ---------------------------------------------------------------------------- %
% Plot the pretty pictures!
% ---------------------------------------------------------------------------- %
%INI.finess
%INI.finess.output_dir

%descriptor = sscanf(INI.finess.output_dir,'%s');
%descriptor = descriptor(7:end)

% n1 = frame number
%fname = strcat( strcat( 'density-schl-dt05_800x200_nf', num2str(n1, '%02d' ) ), '.jpg' );
%print(2, '-djpeg', fname );
%fname = strcat( strcat( 'density-contour-dt05_800x200_nf', num2str(n1, '%02d' )), '.jpg' );
%print(3, '-deps', fname  );

