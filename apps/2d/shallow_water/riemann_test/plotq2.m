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
if( m==1 )
  pcolor(xc,yc,qsoln(:,:,m)+aux(:,:,1));
else
  pcolor(xc,yc,qsoln(:,:,m));
end

shading flat;
%yrbcolormap
axis on; box on; grid off;
axis('equal');
axis([-0.01 1.01 -0.01 1.01]);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:0.25:1);
set(gca,'fontsize',16);

if(m==1)
  t1 = title(['Total height at t = ',num2str(time),'     [DoGPack]']); 
else
  t1 = title(['q(',num2str(m),') at t = ',num2str(time),'     [DoGPack]']); 
end
set(t1,'fontsize',16);

% Show a plot of every row (column) of the solution.  For this Riemann problem,
% each row (column) should be identical.
figure(2);
clf;
if( OPT==1 )
    pz=plot(reshape(xc,mx*my,1), reshape(qsoln(:,:,1)+aux(:,:,1),mx*my,1),'bo');
else
    pz = plot(reshape(yc,mx*my,1), reshape(qsoln(:,:,1)+aux(:,:,1),mx*my,1),'bo');
%   imid = floor( mx/2 );
%   pz = plot(yc(imid,:),  qsoln(imid, :, 1), 'bo');
end
set(pz,'linewidth',2);
set(pz,'markersize',8);
hold on;
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:0.5:3.5);
axis([0 1 0 3.5]);
t1 = title(['Total height at t = ',num2str(time),'     [DoGPack]']);

figure(3);
clf;
if( OPT==1 )
    pz = plot(xc(:,1), qsoln(:,1,1)+aux(:,1,1),'bo');
else
    pz = plot(yc(1,:), qsoln(1,:,1)+aux(1,:,1),'bo');
set(pz,'linewidth',2);
set(pz,'markersize',8);
    hold on;
    pz = plot(yc(2,:), qsoln(2,:,1)+aux(2,:,1),'bo');
set(pz,'linewidth',2);
set(pz,'markersize',8);
    pz = plot(yc(3,:), qsoln(3,:,1)+aux(2,:,1),'bo');
    hold off;
end
set(pz,'linewidth',2);
set(pz,'markersize',8);
hold on;
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:0.5:3.5);
axis([0 1 0 3.5]);
t1 = title(['Total height at t = ',num2str(time),'     [DoGPack]']);
    
figure(1)
