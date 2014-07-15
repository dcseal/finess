function fp = buckley_fp( q )

    M = 1/3;
    fp = (2.*M*q.*(1.-q) ) ./ (q.^2 + M*(1.-q).^2).^2;

end


