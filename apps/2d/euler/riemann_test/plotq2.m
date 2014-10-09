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

% gas constant
global INI;
gamma_gas = sscanf(INI.euler.gamma, '%e');

figure(1);
clf;
pcolor(xc,yc,qsoln(:,:,m));
shading flat;
yrbcolormap
axis on; box on; grid off;
axis('equal');
axis([-0.01 1.01 -0.01 1.01]);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:0.25:1);
set(gca,'fontsize',16);
t1 = title(['q(',num2str(m),') at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);
    
fids  = fopen(['../../../1d/euler/shock_tube/output/qhelp.dat'],'r');
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

  fname = ['../../../1d/euler/shock_tube/output/',...
         num2str(n1+10000),'.dat'];
  fname(37) = 'q';
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
  if (OPT==1)
    pz=plot(reshape(xc,mx*my,1),reshape(qsoln(:,:,1),mx*my,1),'bo');
  else
    pz=plot(reshape(yc,mx*my,1),reshape(qsoln(:,:,1),mx*my,1),'bo');
  end
  set(pz,'linewidth',2);
  set(pz,'markersize',8);
  hold on;
  pr = plot(xc_1d,qsoln_1d(:,1),'r-');
  set(pr,'linewidth',2);
  hold off;
  axis on; box on; grid off;
  set(gca,'plotboxaspectratio',[1.5 1 1]);
  set(gca,'fontsize',16);
  set(gca,'xtick',0:0.25:1);
  set(gca,'ytick',0:0.5:3.5);
  axis([0 1 0 3.5]);
  t1 = title(['Density at t = ',num2str(time),'     [DoGPack]']);
end

figure(1)
