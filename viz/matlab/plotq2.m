%PLOTQ2    User supplied plotting script.
%
% The script PLOTQ2 is a user supplied plotting script, with copies of
% *this* script formatted and placed in each application.  
%
% The file you are currently reading (located in $FINESS/viz/matlab/plotq2.m)
% is the default script, and if an application doesn't define a replacement
% for this script, this is the one that will be run.
%
% Basic information you have available to you:
%
%     Basic information:
%                  mx, my:  number of points in each coordinate direction
% [xlow,xhigh,ylow,yhigh]:  min/max values of grid
%                    meqn:  number of equations
%                    maux:  number of aux components
%
%   Grid information:
%       (xc,yc): grid points (cell centers), size = (mx,my)
%
%   Solution information:
%         qsoln:  solution sampled on mesh, size = (mx,my,meqn)
%           aux:  aux components sampled on mesh, size = (mx,my,maux)
%
% Note: Be careful about redefining variables that are used in plotq2.m.
% For example, q is used as a counter, so do not set something like "q = 5;"
%
% See also: PLOTFIN1.

figure(1);
clf;
surf(xc, yc, qsoln(:, :, m));
yrbcolormap
hold on;
plot([xlow xhigh xhigh xlow xlow],[ylow ylow yhigh yhigh ylow],'k-');
hold off;
axis on; box on; grid off;
axis('equal');
axis([xlow-xeps xhigh+xeps ylow-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['q(',num2str(m),') at t = ',num2str(time),'     [FINESS]']); 
set(t1,'fontsize',16);
c1=colorbar;
set(c1,'fontsize',16);
shading flat;

figure(2);
clf;
contour(xc,yc,qsoln(:,:,m),15,'k');
hold on;
plot([xlow xhigh xhigh xlow xlow],[ylow ylow yhigh yhigh ylow],'k-');
hold off;
axis on; box on; grid off;
axis('equal');
axis([xlow-xeps xhigh+xeps ylow-yeps yhigh+yeps]);
set(gca,'fontsize',16);
t1 = title(['q(',num2str(m),') at t = ',num2str(time),'     [FINESS]']); 
set(t1,'fontsize',16);

figure(1);
