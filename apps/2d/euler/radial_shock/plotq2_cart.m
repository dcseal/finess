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
if( exist( [outputdir,'/eulerhelp.dat'], 'file' ) )
    fids  = fopen([outputdir,'/eulerhelp.dat'],'r');
    if(fids==-1)
      error(['File  ',outputdir,'/eulerhelp.dat  not found.']);
    end
    gamma_gas  = fscanf(fids,'%e',1);
    fclose(fids);
else
    gamma_gas = 1.4;
    disp(['Setting gamma_gas = ', num2str( gamma_gas, '%f' ) ]);
end

rr = reshape(sqrt( xc.^2 + yc.^2 ),mx*my,1);
    
rho_r = reshape(qsoln(:,:,1),mx*my,1);
ux_r  = reshape(qsoln(:,:,2),mx*my,1)./rho_r;
uy_r  = reshape(qsoln(:,:,3),mx*my,1)./rho_r;
uz_r  = reshape(qsoln(:,:,4),mx*my,1)./rho_r;
E_r   = reshape(qsoln(:,:,5),mx*my,1);
u_r   = (reshape(xc,mx*my,1).*ux_r + reshape(yc,mx*my,1).*uy_r)./rr;
p_r   = (gamma_gas-1).*(E_r - 0.5.*rho_r.*(ux_r.^2 + uy_r.^2 + uz_r.^2));

figure(1);
clf;
pcolor(xl,yl,qaug(:,:,m));
shading flat;
yrbcolormap
axis on; box on; grid off;
axis('equal');
axis([-1.05 1.05 -1.05 1.05]);
set(gca,'xtick',-1:0.5:1);
set(gca,'ytick',-1:0.5:1);
set(gca,'fontsize',16);
t1 = title(['q(',num2str(m),') at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
    

fids  = fopen(['../../../1d/euler/shock_tube_radial/output/qhelp.dat'],'r');
if (fids>0)
  nplot_1d = fscanf(fids,'%d',1);
  meqn_1d  = fscanf(fids,'%d',1);
  maux_1d  = fscanf(fids,'%d',1);
  meth1_1d = fscanf(fids,'%d',1);
  mx_1d    = fscanf(fids,'%d',1);
  xlow_1d  = fscanf(fids,'%e',1);
  xhigh_1d = fscanf(fids,'%e',1);
  dx_1d    = fscanf(fids,'%e',1);
  fclose(fids);
  
  mx_old_1d = mx_1d;
  mx_1d = mx_1d*points_per_dir;
  dx_1d = (xhigh_1d-xlow_1d)/mx_1d;
  xc_1d = transpose(linspace(xlow_1d+dx_1d/2,xhigh_1d-dx_1d/2,mx_1d));

  phi_1d = SampleBasis1(points_per_dir,meth1_1d);

  fname = ['../../../1d/euler/shock_tube_radial/output/',...
           num2str(n1+10000),'.dat'];
  fname(44) = 'q';
  fids = fopen(fname,'r');
  time_1d = fscanf(fids,'%e',1);
  qtmp = fscanf(fids,'%e',[1,inf]);
  fclose(fids);
  qtmp = transpose(qtmp);
  qtmp  = reshape(qtmp,mx_old_1d,meqn_1d,meth1_1d);
  qsoln_1d = zeros(mx_1d,meqn_1d);
  for i=1:mx_old_1d
    for me=1:meqn_1d
      for ii=1:points_per_dir
        v1(1:meth1_1d,1) = phi_1d(ii,:);
        v2(1:meth1_1d,1) = qtmp(i,me,:);
        qsoln_1d((i-1)*points_per_dir+ii,me) = transpose(v1)*v2;
      end
    end
  end
  clear qtmp;
        
  figure(2);
  clf;
  pz=plot(rr,rho_r,'bo');
  set(pz,'linewidth',1);
  set(pz,'markersize',8);
  hold on;
  pr=plot(xc_1d,qsoln_1d(:,1),'r-');
  set(pr,'linewidth',2);
  hold off;
  axis on; box on; grid off;
  axis([0 1 0.5 6.0]);
  set(gca,'plotboxaspectratio',[1.5 1 1]);
  set(gca,'xtick',0:0.25:1);
  set(gca,'ytick',0:1:6);
  set(gca,'fontsize',16);
  t1 = title(['Density at t = ',num2str(time),'     [DoGPack]']); 
  set(t1,'fontsize',16);
  
  figure(3);
  clf;
  press_1d = (gamma_gas-1).*(qsoln_1d(:,5)-0.5*(qsoln_1d(:,2).^2 + ...
                                                qsoln_1d(:,3).^2 + qsoln_1d(:,4).^2)./qsoln_1d(:,1));
  
  pz=plot(rr,p_r,'bo');
  set(pz,'linewidth',1);
  set(pz,'markersize',8);
  hold on;
  pr=plot(xc_1d,press_1d,'r-');
  set(pr,'linewidth',2);
  hold off;
  axis on; box on; grid off;
  axis([0 1 0.6 1.8]);
  set(gca,'plotboxaspectratio',[1.5 1 1]);
  set(gca,'xtick',0:0.25:1);
  set(gca,'ytick',0.6:0.3:1.8);
  set(gca,'fontsize',16);
  t1 = title(['Pressure at t = ',num2str(time),'     [DoGPack]']); 
  set(t1,'fontsize',16);
  
  figure(4);
  clf;
  pz=plot(rr,u_r,'bo');
  set(pz,'markersize',8);
  set(pz,'linewidth',1);
  hold on;
  pr=plot(xc_1d,qsoln_1d(:,2)./qsoln_1d(:,1),'r-');
  set(pr,'linewidth',2);
  hold off;
  axis on; box on; grid off;
  axis([0 1 -0.1 0.25]);
  set(gca,'plotboxaspectratio',[1.5 1 1]);
  set(gca,'xtick',0:0.25:1);
  set(gca,'ytick',-0.1:0.1:0.3);
  set(gca,'fontsize',16);
  t1 = title(['u^r(x,t) at t = ',num2str(time),'     [DoGPack]']); 
  set(t1,'fontsize',16); 
end

figure(1)
