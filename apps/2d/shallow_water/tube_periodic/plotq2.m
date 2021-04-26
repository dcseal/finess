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

OPT = 1;   % if OPT==1, shock aligned in x-direction
           % if OPT==2, shock aligned in y-direction

figure(1);
clf;
pcolor(xc,yc,qsoln(:,:,m));

shading flat;
%yrbcolormap
colormap parula
colorbar
caxis([0.1 1]);
axis on; box on; grid off;
axis('equal');
axis([-0.01 1.01 -0.01 1.01]);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:0.25:1);
set(gca,'fontsize',16);

t1 = title(['Total height at t = ',num2str(time),'     [FINESS]']); 
set(t1,'fontsize',16);

figure(2);
clf;
jmid = floor(my/2);
pz = plot(xc(:,jmid), qsoln(:,jmid,1),'bo');
set(pz,'linewidth',2);
set(pz,'markersize',8);
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:0.5:3.5);
axis([0 1 0 1.1]);
t1 = title(['Total height at t = ',num2str(time),'     [FINESS]']);
    
figure(1)
