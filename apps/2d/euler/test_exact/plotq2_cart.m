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
  fids  = fopen([outputdir,'/eulerhelp.dat'],'r');
  if fids==-1
    error(['File  ',outputdir,'/eulerhelp.dat  not found.']);
  end
  gamma_gas  = fscanf(fids,'%e',1);
  fclose(fids);
  
  % Exact solution:
  %    rho = 1 + 0.2*sin( pi*(x+y) - pi*(u+v)*time )
  %
  % with initial conditions u = 0.7 and v = 0.3, u+v=1, and therefore, the exact
  % solution reduces to:
  rho_ex = 1.0 + 0.2*sin( pi*(xc+yc) - pi*time );
  
  figure(1);
  clf;
  pcolor(xl,yl,qaug(:,:,m));
  shading flat;
  yrbcolormap
  axis on; box on; grid off;
  axis('equal');
  axis([-0.05 2.05 -0.05 2.05]);
  set(gca,'xtick',-4:0.5:4);
  set(gca,'ytick',-4:0.5:4);
  set(gca,'fontsize',16);
  t1 = title(['\rho(t = ', num2str(time,'%2.2e'), ', x)     [DoGPack]']); 
  set(t1,'fontsize',16);
  colorbar;
  caxis([0.5,1.5]);
  
  % Take a slice of the solution at x = 1.0:
  if mod(mx,2) == 0 % even number of cells printed!
    xcenter_index = mx / 2;
  else % an odd number of cells
    xcenter_index = (mx+1) / 2;
  end
  qslice = qsoln(xcenter_index, :, 1)';
  rhoslice_ex = 1.0 + 0.2*sin( pi*(xc(xcenter_index,1)+yc(xcenter_index,:)) - pi*time );

  figure(2)
  clf;
  plot(yc(xcenter_index,:), qslice, 'bo' );
  hold on;
  plot(yc(xcenter_index,:), rhoslice_ex, '-k' );
  
  t1 = title(['q(1) at x = ', num2str(xc(xcenter_index,5), '%1.3f'), ' and t = ',num2str(time),'     [DoGPack]']); 

  set(gca,'fontsize',16);
  set(t1,'fontsize',16);

  err = norm( qsoln(:,:,1) - rho_ex, 2 ) / norm( rho_ex, 2 );
  disp([['L2-error = ', num2str( err, '%2.5e' )]]);
  
  % n1 = frame number
  %fname = strcat( strcat( 'density', num2str(n1, '%02d' ) ), '.jpg' );
  %print(1, '-djpeg', fname );
  %print(3, '-dpdf', 'euler-density-contour.pdf' );
  
