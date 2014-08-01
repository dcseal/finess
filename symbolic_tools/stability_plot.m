% plotting parameters
nu0 = linspace(-0.01, 1.5);
k0  = linspace( -800, 800, 1000 );

[nu,k] = meshgrid( nu0, k0 );

% Constants
I  = sqrt( -1. );
dx = 1.0;

s = nu.^3.*(-exp(4*I*dx*k)/1440 + 13*exp(3*I*dx*k)/720 - 55*exp(2*I*dx*k)/432 + 71*exp(I*dx*k)/540 + 91/360 - 583*exp(-I*dx*k)/1080 + 731*exp(-2*I*dx*k)/2160 - exp(-3*I*dx*k)/12 + 47*exp(-4*I*dx*k)/4320 - exp(-5*I*dx*k)/2160) + nu.^2.*(exp(4*I*dx*k)/480 - 3*exp(3*I*dx*k)/80 + 11*exp(2*I*dx*k)/72 + 61*exp(I*dx*k)/360 - 41/80 - exp(-I*dx*k)/180 + 121*exp(-2*I*dx*k)/360 - exp(-3*I*dx*k)/8 + 31*exp(-4*I*dx*k)/1440 - exp(-5*I*dx*k)/720) + nu.*(exp(2*I*dx*k)/20 - exp(I*dx*k)/2 - 1/3 + exp(-I*dx*k) - exp(-2*I*dx*k)/4 + exp(-3*I*dx*k)/30) + 1;

sfe = (1.0-nu) + nu.*exp(-I*k*dx);

figure(1);
clf;

[C,h] = contour( nu, k, abs( s ), [1,1] );

xlabel( 'CFL number', 'FontSize', 16 );
ylabel( 'Wave number', 'FontSize', 16 );

set(gca,'xtick', -2:0.05:2);
%set(gca,'ytick',-2:0.25:2);

t1 = title('Von Neumann stability analysis' );
set(t1, 'fontsize', 16 );

