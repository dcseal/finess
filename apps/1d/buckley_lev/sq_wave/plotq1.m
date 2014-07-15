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
pz=plot(xc,qsoln(:,m),'bo');
hold off;
axis on; box on; grid off;
axis([-1 1 -0.1 1.1+0.1]);
set(gca,'plotboxaspectratio',[2 1 1]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.5:2);
set(gca,'fontsize',16);
grid on;
t1 = title(['q(x,t) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

%% -- Add in the exact solution -- %%
qs_left_problem  = 0.1339745962155613;
qs_right_problem = 0.5;

qsample = linspace(qs_right_problem, 1.0 );
dt      = time;

ql = linspace(0., qs_left_problem );
left_shock = -0.5 + dt * buckley_fp( qs_left_problem );
fl = -0.5 + dt*buckley_fp( ql );

% segment 1
hold on;
plot( [-1 -0.5], [0 0], '-r' )
hold on;

% segment 2
plot( fl, ql, '-r' );
hold on;

% segment 3
plot( fl, ql, '-r' );
plot( [left_shock left_shock], [qs_left_problem 1.0], '-r' );
hold on;

% segment 4
plot( [left_shock 0], [1 1], '-r' );


% segment 5
plot( dt*buckley_fp(qsample), qsample, '-r'  );
hold on;

% segment 6
plot( [ dt*buckley_fp(qs_right_problem) dt*buckley_fp(qs_right_problem) ], ...
      [ qs_right_problem 0 ], '-r' );
hold on;

% segment 7
plot( [dt*buckley_fp(qs_right_problem) 1.0], [0 0], '-r' );

hold on;
plot( [-1 1], [0 0], '--k' );

%axis( [-1 1 -0.1 1.1] );
hold off;


figure(2);
clf;
pz=plot(xc,qsoln(:,m),'bo');
hold off;
axis on; box on; grid off;
axis([-0.1 0.1 0.9 1.1]);
set(gca,'plotboxaspectratio',[2 1 1]);
set(gca,'xtick',-2:0.05:2);
set(gca,'ytick',-2:0.05:2);
grid on;
set(gca,'fontsize',16);
t1 = title(['q(x,t) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);

%% -- Add in the exact solution -- %%
qs_left_problem  = 0.1339745962155613;
qs_right_problem = 0.5;

qsample = linspace(qs_right_problem, 1.0 );
dt      = time;

ql = linspace(0., qs_left_problem );
left_shock = -0.5 + dt * buckley_fp( qs_left_problem );
fl = -0.5 + dt*buckley_fp( ql );

% segment 1
hold on;
plot( [-1 -0.5], [0 0], '-r' )
hold on;

% segment 2
plot( fl, ql, '-r' );
hold on;

% segment 3
plot( fl, ql, '-r' );
plot( [left_shock left_shock], [qs_left_problem 1.0], '-r' );
hold on;

% segment 4
plot( [left_shock 0], [1 1], '-r' );


% segment 5
plot( dt*buckley_fp(qsample), qsample, '-r'  );
hold on;

% segment 6
plot( [ dt*buckley_fp(qs_right_problem) dt*buckley_fp(qs_right_problem) ], ...
      [ qs_right_problem 0 ], '-r' );
hold on;

% segment 7
plot( [dt*buckley_fp(qs_right_problem) 1.0], [0 0], '-r' );

hold on;
plot( [-1 1], [0 0], '--k' );

%axis( [-1 1 -0.1 1.1] );
hold off;
