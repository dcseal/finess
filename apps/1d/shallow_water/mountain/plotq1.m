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

figure(1);
clf;
pz=plot(xc,qsoln(:,1),'bo');
set(pz,'markersize',8);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
t1 = title(['Height at t = ',num2str(time),'     [FINESS]']);
axis([0 1 0.25 1.15]);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',-1:0.1:2);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(t1,'fontsize',16);

figure(2);
clf;  
pz=plot(xc, qsoln(:,2)./qsoln(:,1), 'bo');
set(pz,'markersize',8);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
t1 = title(['Velocity at t = ',num2str(time),'     [FINESS]']);  
axis([0 1 -0.1 0.1]);
axis auto;
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',-1:0.1:1.5);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(t1,'fontsize',16);



