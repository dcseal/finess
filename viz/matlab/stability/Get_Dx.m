function [Dx,Dxx]  = Get_Dx( mx )
%Get matrices for centered finite difference approximations.
%
% This routine returns a (sparse) matrix using the stencil of minimal size.
%
% See also: Get_Qt.

    dx = 1./mx;

    % Stencil (with linear weights)
    Lx  = [1 -8 0 8 -1]/(12*dx);
%   Lx  = [-2 -1 10 1 2];

    Lxx = [-1 16 -30 16 -1]/(12*dx*dx);
%   Lxx = [-2 -1 10 1 2];

    % Set up the matrix
    e  = ones( mx, 1 );
    Dx = spdiags( e*Lx, [-2 -1 0 1 2], mx, mx );

    % Periodic boundary conditions
    Dx(1, mx-1:mx ) = Lx(1:2);
    Dx(2,      mx ) = Lx(1  );
    Dx(mx-1, 1   )  = Lx(5);
    Dx(mx,   1:2 )  = Lx(4:5);

    % Second derivative 
    Dxx = spdiags( e*Lxx, [-2 -1 0 1 2], mx, mx );

    % Periodic boundary conditions
    Dxx(1, mx-1:mx ) = Lxx(1:2);
    Dxx(2,      mx ) = Lxx(1  );
    Dxx(mx-1, 1   )  = Lxx(5);
    Dxx(mx,   1:2 )  = Lxx(4:5);

end
