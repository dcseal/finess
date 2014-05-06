function A = Get_Qt( mx )
% Compute the time derivative of the large system.
%
% See "On the Linear Stability of the Fifth-Order WENO Discretization", Journal
% of Scientific Computing May 2011, Volume 47, Issue 2, pp 127-149.

    dx = 1./mx;

    % Stencil (with linear weights)
%   fiph  = [2/60 -13/60 47/60 27/60 -3/60];
%   Lp    = [0 fiph];
%   Lm    = [fiph 0];
%   L     = -(Lp-Lm)/dx;
    L = -[-2/60 15/60 -1 20/60 30/60 -3/60]/dx;

    % Set up the matrix
    e  = ones( mx, 1 );
    A = spdiags( e*L, [-3 -2 -1 0 1 2], mx, mx );

    % Periodic boundary conditions
    A(1, mx-2:mx ) = L(1:3);
    A(2, mx-1:mx ) = L(1:2);
    A(3,      mx ) = L(1  );
    A(mx-1, 1   )  = L(6);
    A(mx,   1:2 )  = L(5:6);

end
