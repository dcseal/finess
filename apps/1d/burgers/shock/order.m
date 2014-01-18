clear;
fnt = 16;
format short e;

fid = fopen('err2.dat','r');
err_vec = fscanf(fid,'%f %f',[2,inf]);
fclose(fid);
err_vec = err_vec';
h1 = err_vec(:,1);
E1 = err_vec(:,2);

fid = fopen('err3.dat','r');
err_vec = fscanf(fid,'%f %f',[2,inf]);
fclose(fid);
err_vec = err_vec';
h2 = err_vec(:,1);
E2 = err_vec(:,2);

% Least squares
tmp = ones(size(h1));
A = [tmp, log(h1)];
B = inv(((A')*A));
tmp1 = B*(A')*log(E1);
c1 = exp(tmp1(1,1));
p1 = tmp1(2,1);
pmod1 = round(100*p1)/100;

% Least squares
tmp = ones(size(h2));
A = [tmp, log(h2)];
B = inv(((A')*A));
tmp2 = B*(A')*log(E2);
c2 = exp(tmp2(1,1));
p2 = tmp2(2,1);
pmod2 = round(100*p2)/100;

% Plot
figure(1)
clf

ttm2=loglog(h1,c1*h1.^(p1),'r-');
set(ttm2,'linewidth',2);

hold on;

ttm2=loglog(h2,c2*h2.^(p2),'r-');
set(ttm2,'linewidth',2);

ttm1=loglog(h1,E1,'bo');
set(ttm1,'markersize',8);
set(ttm1,'linewidth',4);

ttm1=loglog(h2,E2,'bo');
set(ttm1,'markersize',8);
set(ttm1,'linewidth',4)

t1=title('1-norm error');
set(t1,'fontsize',fnt);
x1=xlabel('step size');
y1=ylabel('error');
set(x1,'fontsize',fnt);
set(y1,'fontsize',fnt);
set(gca,'fontsize',fnt);
axis([5e-4 5e-2 1e-6 1e+1]);
set(gca,'plotboxaspectratio',[1.5 1 1])
set(gca,'xtick',[1e-3,1e-2,1e-1]);
set(gca,'ytick',[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e-0,1e+1,1e+2]);

ta1=gtext([ 'p = ',num2str(pmod1)]);
set(ta1,'fontsize',fnt+2);

ta1=gtext([ 'p = ',num2str(pmod2)]);
set(ta1,'fontsize',fnt+2);

hold off