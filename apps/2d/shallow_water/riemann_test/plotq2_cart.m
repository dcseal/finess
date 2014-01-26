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

OPT = 2;   % if OPT==1, shock aligned in x-direction
           % if OPT==2, shock aligned in y-direction

figure(1);
clf;
if( m==1 )
  pcolor(xl,yl,qaug(:,:,m)+aux_aug(:,:,1));
else
  pcolor(xl,yl,qaug(:,:,m));
end

shading flat;
yrbcolormap
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
    
figure(1)
