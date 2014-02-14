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
%set(pz,'linewidth',2);
hold off;
axis on; box on; grid off;
%axis([0 1 -0.6 1.6]);
axis([0 2 -1.1+0.5 1.1+0.5]);
set(gca,'plotboxaspectratio',[2 1 1]);
set(gca,'xtick',-2:0.25:2);
set(gca,'ytick',-2:0.5:2);
set(gca,'fontsize',16);
t1 = title(['q(x,t) at t = ',num2str(time),'     [DoGPack]']); 
set(t1,'fontsize',16);


%sp = zeros(mx_old,1);
%se = zeros(mx_old,1);
%flag = zeros(mx_old,1);
%for i=1:mx_old
%  Q2sum = 0.0;
%      
%  for k=1:(meth1-2)
%    Q2sum = Q2sum + qcoeffs(i,1,k)^2;
%  end
%      
%  Q2kmax = qcoeffs(i,1,meth1-1)^2;
%  Q2sum = Q2sum + Q2kmax;
%
%  if Q2kmax>1.0e-10
%    sp(i) = (Q2kmax/Q2sum);
%  else
%    sp(i) = 0.0;
%  end
%  
%  Q2kmax = qcoeffs(i,1,meth1)^2;
%  Q2sum = Q2sum + Q2kmax;
%  
%  if Q2kmax>1.0e-10
%    se(i) = (Q2kmax/Q2sum);
%  else
%    se(i) = 0.0;
%  end
%  
%  if (se(i)>=(0.1*(meth1)^(-4)))
%    flag(i) = 1;
%  end
%  
%end

%figure(2);
%clf;
%pz=plot(xc_old,se,'bo');
%%hold on
%%pzz=plot(xc_old,sp,'rx');
%hold off;
%axis on; box on; grid off;
%set(gca,'plotboxaspectratio',[2 1 1]);
%set(gca,'fontsize',16);

%figure(3);
%clf;
%pz=plot(xc_old,flag,'bo');
%hold off;
%axis on; box on; grid off;
%set(gca,'plotboxaspectratio',[2 1 1]);
%set(gca,'fontsize',16);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for the exact solution.
%
% You can comment out this entire section of the script if you're not
% computing an exact solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Parameters used for computing the exact solution
 t = time;
 q0 = @(z)(0.5 + sin(pi*z));
 qp = @(z)(pi*cos(pi*z));
 f  = @(z)(q0(z)*t + z - xc);
 fp = @(z)(qp(z)*t + 1.0);

 % Convergence parameters:
 tol = 1e-14;
 max_iter = 1000;
 tfinal = 0.5/pi;

 % Only check for exact solution at the final time
 if( abs(t-tfinal) < 1e-10 )
%  disp('  Computing exact solution: comment out this section of plotq1 if not desired');
   xi  = xc;   % initial guess
   err = max(abs( f(xi) ) );

   % Newton iteration til convergence:
   num_iter = 0;
   while( err > tol )
   	xi = xi - (f(xi) ./ fp(xi));
   	err = max(abs(f(xi)));
 	num_iter = num_iter+1;
 	if( num_iter > max_iter )
 		disp('****   too many iterations ****' );
        disp('you are computing junk');
 		break;
 	end
   end

   qex = q0(xc-t*q0(xi));
%  disp(['   max(qex-q0(xi)) = ', num2str(max(qex-q0(xi)),'%0.8e')]);
 
   hold on;
   pt = plot(xc,qex, 'r-');
   set(pt,'linewidth',1);

%  disp(['   min(xi) = ', num2str(min(xi)), '    max(xi) = ', num2str(max(xi))]);

    % print error in exact vs computed solution
    err1 = norm(qsoln-qex,1)/norm(qex,1);
    err2 = norm(qsoln-qex,2)/norm(qex,2);
    erri = max( abs( qsoln-qex ) ) / max( abs( qex ) );

    % Friendly helper message: (commented to pull convergence numbers more easily)
    fprintf(1, '%d %2.4e %2.4e %2.4e\n', mx, err1, err2, erri );

    fid = fopen('errors.dat', 'a' );
    fprintf(fid, '%d %2.15e %2.15e %2.15e\n', mx, err1, err2, erri );
    fclose( fid );

    % Plot the error
%   figure(4);
%   clf;
%   pz=plot(xc, qsoln(:,m)-qex, 'bo');
%   set(pz,'linewidth',2);
%   hold off;
%   axis on; box on; grid off;
%   axis auto;
%   set(gca,'plotboxaspectratio',[2 1 1]);
%   set(gca,'xtick',-2:0.25:2);
%   set(gca,'ytick',-2:0.5:2);
%   set(gca,'fontsize',16);
%   t1 = title(['q(x,t) - qex(x,t) at t = ',num2str(time),'     [DoGPack]']); 
%   set(t1,'fontsize',16);

 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
