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

h = qsoln(:,:,1);
u = qsoln(:,:,2)./h;
v = qsoln(:,:,3)./h;

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

t1 = title(['Total height at t = ',num2str(time)]);
%t1 = title(['Total height at t = ',num2str(time),'     [FINESS]']); 
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
grid on;
t1 = title(['Total height at t = ',num2str(time),'     [FINESS]']);
t1 = title(['Total height at t = ',num2str(time), ', y = ', num2str(yc(1,jmid))]);


% ---------------------------------------------------------------------------- %
% Schlieren plot for water height
% ---------------------------------------------------------------------------- %

figure(3);
clf;

%colors=linspace(1.0,3.0,10);
dx=xc(2,1)-xc(1,1);
[hx,hy]=gradient(h,mx,mx);
vel=sqrt(hx.^2+hy.^2)';

colormap( flipud(gray(2048)).^10 )
hp = pcolor(xc,yc,vel');
set(hp, 'EdgeColor', 'none');

% Fancier labels for the plots
axis on; box on; grid off;
axis('equal');
axis([-0.05 1.05 -0.05 1.05]);
set(gca,'xtick',-4:0.5:4);
set(gca,'ytick',-4:0.5:4);
set(gca,'fontsize',16);
t1 = title(['Height at t = ',num2str(time)]); 
set(t1,'fontsize',16);

figure(1)

% Print the pretty pictures!
if( exist( [outputdir,'/photos'], 'dir' ) )
  folder_name = strcat( outputdir, '/photos/' );
  fname = strcat( strcat('height.mx.', num2str(mx,'%03d'), '.my.', num2str(my,'%03d'), '.nf.', num2str(n1, '%03d') ), '.eps' );
  print(1,'-depsc', strcat(folder_name, fname ) );
  
  fname = strcat( strcat('height-cross-section.mx.', num2str(mx,'%03d'), '.my.', num2str(my,'%03d'), '.nf.', num2str(n1, '%03d') ), '.eps' );  
  print(2,'-depsc', strcat(folder_name, fname ) );

  fname = strcat( strcat('height-schlieren.mx.', num2str(mx,'%03d'), '.my.', num2str(my,'%03d'), '.nf.', num2str(n1, '%03d') ), '.eps' );
  print(3,'-depsc', strcat(folder_name, fname ) );
  
end