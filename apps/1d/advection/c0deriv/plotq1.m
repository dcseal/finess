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

if time==0
  qex = qsoln;
end

figure(1);
clf;
pt=plot(xc,qex(:,1),'r-');
set(pt,'linewidth',1.5);
hold on;
pz=plot(xc,qsoln(:,1),'bo');
set(pz,'markersize',8)
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
axis([-1 1 -0.1 1.1]);
set(gca,'plotboxaspectratio',[2 1 1]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.5:2);
set(gca,'fontsize',16);
t1 = title(['q(x,t) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

if ( abs(time-round(time)) < 1.0e-10 && time > 0)
  err1 = norm(qsoln-qex,1)/norm(qex,1);
  err2 = norm(qsoln-qex,2)/norm(qex,2);
  disp(' ');
  disp(['   dx = ',num2str(dx,'%0.8e'),['         err1 = ' ...
                          ''],num2str(err1,'%0.8e')]);
  disp(' ');
  disp(['   dx = ',num2str(dx,'%0.8e'),['         err2 = ' ...
                      ''],num2str(err2,'%0.8e')]);
  disp(' ');
end

disp([['Min value of q = ', num2str( min(qsoln(:,1)), '%1.10e' ) ]]);
