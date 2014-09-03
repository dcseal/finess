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
%%       (xl,yl): grid points (lower left cell corners), size = (mx+1,my+1)
%%
%%   Solution information:
%%         qsoln:  solution sampled on mesh, size = (mx,my,meqn)
%%           aux:  aux components sampled on mesh, size = (mx,my,maux)
%%          qaug:  solution sampled on mesh, with zero padding to
%%                 make it compatible with surf and pcolor matlab
%%                 plotting tools, size = (mx+1,my+1,meqn)
%%       aux_aug:  aux components sampled on mesh, with zero padding to
%%                 make it compatible with surf and pcolor matlab
%%                 plotting tools, size = (mx+1,my+1,maux)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gas constant
% fids  = fopen([outputdir,'/eulerhelp.dat'],'r');
% if fids==-1
%   disp(['File  ',outputdir,'/eulerhelp.dat  not found.']);
%   disp('Setting gamma = 1.4');
%   gamma_gas = 1.4;
% else
%   gamma_gas  = fscanf(fids,'%e',1);
%   fclose(fids);
% end

%%HACK
gamma_gas = 5/3
%%----

figure(1);
clf;
pcolor(xl,yl,qaug(:,:,m));
shading flat;
yrbcolormap
axis on; box on; grid off;
axis('equal');
axis([-1.05 1.05 -1.05 1.05]);
set(gca,'xtick',-2:1.0:2);
set(gca,'ytick',-2:1.0:2);
set(gca,'fontsize',16);
t1 = title(['t = ',num2str(time),'     [FINESS]']); 
set(t1,'fontsize',16);
colorbar;
caxis auto;

figure(2);
clf;
surf(xl,yl,qaug(:,:,m));
shading flat;
yrbcolormap
axis on; box on; grid off;
axis('equal');
axis([-1.05 1.05 -1.05 1.05]);
set(gca,'xtick',-2:1.0:2);
set(gca,'ytick',-2:1.0:2);
set(gca,'fontsize',16);
t1 = title(['t = ',num2str(time),'     [FINESS]']); 
set(t1,'fontsize',16);
colorbar;
caxis auto;



% Grey scale figure plot that you see in a lot of papers.
figure(3);
clf;
contour(xl, yl, qaug(:,:,m), linspace(-2, 1.2, 20), '-k' );
axis on; box on; grid off;
axis('equal');
axis([-1.05 1.05 -1.05 1.05]);
set(gca,'xtick',-2:1.0:2);
set(gca,'ytick',-2:1.0:2);
set(gca,'fontsize',16);
t1 = title(['t = ',num2str(time),'     [FINESS]']); 
set(t1,'fontsize',16);

figure(1)

% n1 = frame number
%fname = strcat( strcat( 'density', num2str(n1, '%02d' ) ), '.jpg' );
%print(1, '-djpeg', fname );
%fname = strcat( strcat( 'density-contour', num2str(n1, '%02d' ) ), '.eps' );
%print(3, '-deps', fname  );

