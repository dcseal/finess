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
gamma = 5/3;
%% ----

figure(1);
clf;
%pz=plot(xc,qsoln(:,m),'r-');
pz=plot(xc,qsoln(:,m),'bo');
set(pz,'linewidth',2);
set(pz,'markersize',8);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
set(gca,'xtick',-1:0.5:1);
if m == 1      
  set(gca,'ytick',0:0.2:1);
  axis([-1 1 0.1 1.1]);
  t1 = title(['Density at t = ',num2str(time),'     [FINESS]']);
elseif m==6
  set(gca,'ytick',0.72:0.01:0.78);
  axis([-1 1 0.72 0.78]);
  t1 = title(['B^1(x,t) at t = ',num2str(time),'     [FINESS]']);
elseif m==7
  set(gca,'ytick',-1:0.5:1);
  axis([-1 1 -1.1 1.1]);
  t1 = title(['B^2(x,t) at t = ',num2str(time),'     [FINESS]']);
    else
      t1 = title(['q(',num2str(m),') at t = ',num2str(time),'     [FINESS]']);
end
set(t1,'fontsize',16);
    
press = (gamma-1)*(qsoln(:,5)-0.5*(qsoln(:,2).^2+qsoln(:,3).^2+ ...
                                   qsoln(:,4).^2)./ ...
                   qsoln(:,1)-0.5*(qsoln(:,6).^2+qsoln(:,7).^2+ ...
                                   qsoln(:,8).^2));
    
figure(2)
clf
%pz=plot(xc,press,'r-');
pz=plot(xc,press,'bo');
set(pz,'linewidth',2);
set(pz,'markersize',8);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
t1=title(['Pressure at t = ',num2str(time),'     [FINESS]']);
set(t1,'fontsize',16);
axis([-1 1 0 1.05]);
