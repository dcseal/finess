clear

mx = 300;
dx = 1./mx;
dt = 0.8*dx;

% Construct matrices for FD approximations
I        = speye( mx );
A        = Get_Qt( mx );
[Dx,Dxx] = Get_Dx( mx );

movie = 1;

% Matrices that result from an RK scheme.  (That is, compose the stencil with
% itself)
%   A2 = A*A;
%   A3 = A2*A;
%   A4 = A3*A;
%   A5 = A4*A;
%   A6 = A5*A;
%   A7 = A6*A;
%   A8 = A7*A;

% Fully discrete matrices from compositions with yourself
%   M1 = I  - dt*A;
%   M2 = M1 - dt^2/2   * A2;
%   M3 = M2 - dt^3/6   * A3;
%   M4 = M3 - dt^4/24  * A4;
%   M5 = M4 - dt^5/120 * A5;
%   M6 = M5 - dt^6/720 * A6;
%   M7 = M6 - dt^7/factorial(7) * A7;
%   M8 = M7 - dt^8/factorial(8) * A8;

%   1.0 - max( abs( eig( full( M1 ) ) ) );
%   1.0 - max( abs( eig( full( M2 ) ) ) );
%   1.0 - max( abs( eig( full( M3 ) ) ) );
%   1.0 - max( abs( eig( full( M4 ) ) ) );
%   1.0 - max( abs( eig( full( M5 ) ) ) );
%   1.0 - max( abs( eig( full( M6 ) ) ) );
%   1.0 - max( abs( eig( full( M7 ) ) ) );
%   1.0 - max( abs( eig( full( M8 ) ) ) );

figure(1);
clf;
plot( eig( full( Get_Qt( 100  ) ) ), 'ko' );
hold on;
%plot( eig( full( Get_Qt( 1000 ) ) )/10, 'go' );
hold off;

if( movie )

    M = I + dt*A*( I - dt/2*Dx + dt^2/6*Dxx );

    mt    = 100;
    numax = 1.2;
    dt_vec = linspace(0, numax*dx, mt );
    for n=1:length(dt_vec)

        dt = dt_vec(n);
        M = I + dt*A*( I - dt/2*Dx + dt^2/6*Dxx );

        vals = abs( eig( full( M ) ) );

        if( max( vals ) > 1.0 )
            disp(['nu = ', num2str( dt/dx )] );
        end

        figure(1);
        clf
        plot( 1:mx, vals, 'go' );
        axis( [0 mx 0, 1.5] );
        hold on;
        plot([0 mx],[1 1], 'k--');
        hold on;
        plot([0 0], [0 1.5], 'k--');
        title(['nu = ', num2str( dt/dx )] );
        hold off;

        % PIF-WENO (localized stencil)
        figure(2);
        clf;
        plot( eig( full( M ) ), 'ko' );
        axis( [0 1 -1, 1] );
        title(['nu = ', num2str( dt/dx )]);

    end

    dt = 1.1*dx;
    M  = I + dt*A*( I - dt/2*Dx + dt^2/6*Dxx );

end
