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
t1 = title(['Height at t = ',num2str(time),'     [DoGPack]']);
axis([0 1 0.75 3.25]);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:1:3);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(t1,'fontsize',16);

%% Add in the exact solution:

% left values:
ul = 0.;
hl = 3.;

% constant intermediate state:
ustar = 7.448542169801264e-01;
hstar = 1.848576603096757e+00;

% right values from initial conditions:
ur = 0.;
hr = 1.;

% feet of 1-rarefaction
left_foot  = 0.5 + time*( ul - sqrt( hl ) );
right_foot = 0.5 + time*( ustar - sqrt( hstar ) );

% location of shock:
s = (hstar*ustar - hr*ur) / (hstar - hr );
xshock = 0.5 + s*time;

% Region 1:
hold on;
plot( [xc(1), left_foot], [hl, hl], '-r', 'LineWidth', 3 );

% rarefaction:
if( time > 0 )

    x_rarefaction = linspace( left_foot, right_foot );

    A = ul + 2.*sqrt(hl);
    ht = 1./9.*( A - (x_rarefaction-0.5)/time ).^2;

    plot( x_rarefaction, ht, '-r', 'LineWidth', 3 );

else
    plot( [0.5, 0.5], [hl, hr], '-r', 'LineWidth', 3 );
end

% region 2 (intermediate state):
plot( [right_foot, xshock], [hstar, hstar], '-r', 'LineWidth', 3 );
plot( [xshock, xshock], [hstar, hr], '-r', 'LineWidth', 3 );

% region 3
plot( [xshock, xc(end)], [hr, hr], '-r', 'LineWidth', 3 );


figure(2);
clf;  
pz=plot(xc,qsoln(:,2),'bo');
set(pz,'markersize',8);
set(pz,'linewidth',1);
hold off;
axis on; box on; grid off;
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(gca,'fontsize',16);
t1 = title(['hu^1(x,t) at t = ',num2str(time),'     [DoGPack]']);  
axis([0 1 -0.2 1.6]);
set(gca,'xtick',0:0.25:1);
set(gca,'ytick',0:0.5:1.5);
set(gca,'plotboxaspectratio',[1.5 1 1]);
set(t1,'fontsize',16);

figure(3);
clf;

h = qsoln(:,1);
u = qsoln(:,2)./h;

pz=plot(xc, u + 2.*sqrt(h), 'k-' );
hold on;
pz=plot(xc, u - 2.*sqrt(h), 'r-' );
hold on;;

%   set(pz,'markersize',8);
%   set(pz,'linewidth',1);
%   hold off;
%   axis on; box on; grid off;
%   set(gca,'plotboxaspectratio',[1.5 1 1]);
%   set(gca,'fontsize',16);
%   t1 = title(['hu^1(x,t) at t = ',num2str(time),'     [DoGPack]']);  
%   axis([0 1 -0.2 1.6]);
%   set(gca,'xtick',0:0.25:1);
%   set(gca,'ytick',0:0.5:1.5);
%   set(gca,'plotboxaspectratio',[1.5 1 1]);
%   set(t1,'fontsize',16);
