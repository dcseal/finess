%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%   Information you have available to you:
%%
%%     Basic information:
%%                  mx, my:  number of points in each coordinate direction
%% [xlow,xhigh,ylow,yhigh]:  min/max values of grid
%%                    meqn:  number of equations
%%                    maux:  number of aux components
%%                   meth1:  spatial order of accuracy
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
if( isfield( INI.euler, 'gamma' ) )
    gamma = sscanf( INI.euler.gamma, '%f' );
else
    disp('Seeting gamma = 1.4');
    gamma = 1.4;
end

figure(1);
clf;
v=linspace(0,6,40);
contourf(xc, yc, qsoln(:, :, m), v);
%pcolor(xc, yc, qsoln(:, :, m), v);
shading flat;
yrbcolormap
axis on; box on; grid off;
%axis('equal');
%axis([-0.05 3.25 -0.05 1.05]);
%set(gca,'xtick',-4:0.5:4);
%set(gca,'ytick',-4:0.5:4);
%set(gca,'fontsize',16);
t1 = title(['Density at t = ',num2str(time),'     [FINESS]']); 
%set(t1,'fontsize',16);
colorbar;
%caxis([1,25]);

% Radial slice of the density
figure(2)
clf
pz=plot(xc(:,1),qsoln(:,1,1),'bo');
set(pz,'linewidth',1);
set(pz,'markersize',6);
%t1 = title(['Density']);
%t1 = title(['\rho(t, x, y)     [FINESS]']); 
t1 = title(['\rho(t = ', num2str(time,'%2.2e'), ', x)     [FINESS]']); 
set(t1,'fontsize',12);
axis on; box on; grid off;

% Exact solution (Pulled from Sedov's book)
r = [0 .1 .15 .2 .25 .3 .3654 .4222 .4748 .518 .5754 .639 .6894 .7274 .7629 .8094 .8442 .8725 .9096 .9295 .9476 .9644 .9802 .998 1];
den = [ 0 0 0 0.0001 .0003 .0008 .0021 .0044 .0079 .0123 .0208 .0362 .0545 .0718 .0975 .1414 .1892 .2427 .3451 .4234 .5164 .6285 .7653 .9973 1];

rho = 1;
d0  = (1.4+1)/(1.4-1)*rho;
r0  = 1;

dd = den*d0;
rr = r*r0;

rr = [rr r0 1.1];
dd = [dd 1 1 ];

hold on;
plot( rr, dd, '-r' );
hold off;
axis([0 1.1 0 6.5]);

figure(1)
