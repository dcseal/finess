nruns = 10;
for n=1:nruns
    outputdir = ['output', num2str( n-1, '%03d' ) ];
    plotfin1_nostop( 1, outputdir, 1 );     % we think there's only one point ...
end
