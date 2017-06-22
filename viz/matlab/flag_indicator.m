function flag = flag_indicator( q , params )

    h = params.dx;
    n = params.mx;
    flag = zeros( size(q) );

    discont_margin = params.eps_indicator*h;  %1.1*h^(3/2);

    for i = 2:1:n-1
        d = abs( q(i) - 0.5*( q(i-1) + q(i+1) ) );
        if d > discont_margin
            flag(i) = 1;
        else
            flag(i) = 0;
        end 
    end

end
