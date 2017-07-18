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

% close all;    % Close any open figures

% gas constant
gamma_gas = sscanf(INI.euler.gamma, '%e');

% Riemann problem number.  
%
% Here, we follow the naming convention found in 
% "A posteriori subcell limiting of the discontinuous Galerkin finite element method 
% for hyperbolic conservation laws," Dumbser et al, JCP (2014).  The final
% times reported in that paper are the following:
%
% RP[number] | T-final
% --------------------
%        1   | 0.25  (This is shorter than T = 0.8 found in many other papers)
%        2   | 0.25
%        3   | 0.3
%        4   | 0.25
%        5   | 0.25
rpnum     = sscanf(INI.euler.riemann_problem_number, '%d' );

rho=qsoln(:,:,1);

% figure('Position', [100, 100, 1049, 895]);
% clf;
% pcolor(xc,yc,rho);
% axis([-0.5 0.5 -0.5 0.5])
% %contourf(xc, yc, rho,40 );
% min(min(rho))
% shading flat;
% yrbcolormap
% axis on; box on; grid off;
% axis('equal');
% %axis([-0.05 13.25 -0.05 11.05]);
% %set(gca,'xtick',-4:1.0:14);
% %set(gca,'ytick',-4:1.0:14);
% %set(gca,'fontsize',16);
% t1 = title(['Density at t = ',num2str(time)]); 
% set(t1,'fontsize',16);
% colorbar;
% %caxis([0.8,10]);
% caxis auto;
% %export_fig -transparent Density1.png

% Grey scale figure plot that you see in a lot of papers.
% figure('Position', [100, 100, 1049, 895]);
figure(1);
clf;
%colors=linspace(1.0,3.0,10);
h=xc(2,1)-xc(1,1);
[px,py]=gradient(rho,h,h);
v=sqrt(px.^2+py.^2)';

fprintf('Smallest and largest value of the computed gradient = %f %f\n', max(max(v)), min(min(v)) );
colormap(flipud(gray(2048)).^10)
h = pcolor(xc,yc,v');
set(h, 'EdgeColor', 'none');
%contourf(xc, yc, rho,10, '-k' );

%colorbar()
axis on; box on; grid off;
axis('equal');
%axis([-0.05 13.25 -0.05 11.05]);
%set(gca,'xtick',-4:1.0:14);
%set(gca,'ytick',-4:1.0:14);
%set(gca,'fontsize',16);
t1 = title(['Density at t = ',num2str(time)]); 
axis([-0.5 0.5 -0.5 0.5])       % TODO - read in RP number?
set(t1,'fontsize',16);

% TODO - need to find this file
% from somewhere.  Scott, did you write this yourself or yank it from the web?
% -DS
% export_fig -transparent DensityBW.png     

% figure('Position', [100, 100, 1049, 895]);
  figure(2);
  clf;
  rho=qsoln(:,:,1);
  e=qsoln(:,:,5);
  u1=qsoln(:,:,2)./qsoln(:,:,1);
  u2=qsoln(:,:,3)./qsoln(:,:,1);
  u3=qsoln(:,:,4)./qsoln(:,:,1);
  P=0.4*(e-0.5*rho.*(u1.*u1+u2.*u2+u3.*u3));
  
  c=sqrt(1.4*P./rho);
  M=sqrt(u1.*u1+u2.*u2)./c;
  %contour(xc, yc, P, linspace(0.091, 38, 50), '-k' );
  contourf(xc, yc, P, 40, '-k' );
  colorbar()
  %axis on; box on; grid off;
  %axis('equal');
  %axis([-0.05 13.25 -0.05 11.05]);
  
  axis([-0.5 0.5 -0.5 0.5])
  %set(gca,'xtick',-4:1.0:14);
  %set(gca,'ytick',-4:1.0:14);
  set(gca,'fontsize',16);
  t1 = title(['Pressure at t = ',num2str(time)]); 
  set(t1,'fontsize',16);

% Note: n1 = frame number so you can use this to keep track of each picture
% when coupled with something like plotfin2_nostop (or simply plotfin2 if you
% want to click through each image).

%   descriptor = erase(sscanf(INI.finess.output_dir,'%s'),'output');
    
    descriptor = sscanf(INI.finess.output_dir,'%s');
    descriptor = descriptor(7:end);

%
% fname = strcat( strcat( 'density-dt07-128-frame', num2str(n1, '%02d' ) ), '.jpg' );
  fname = [['density-', descriptor, '-frame', num2str(n1, '%02d' ) , '.jpg' ]];
; print(1, '-djpeg', fname );

% fname = strcat( strcat( 'density-contour', num2str(n1, '%02d' ) ), '.eps' );
% print(3, '-deps', fname  );


