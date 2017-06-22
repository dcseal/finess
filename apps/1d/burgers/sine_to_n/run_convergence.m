% Single script to call all the output folders and perform a convergence
% study.
n = 0;
outputdir = ['output_', num2str( n, '%03d' ), '/parameters.ini' ];
outputdir
while( exist( outputdir, 'file') )
    plotfin1_nostop( outputdir );     % we think there's only one point ...
    n = n+1;
    outputdir = ['output_', num2str( n, '%03d' ), '/parameters.ini' ];
end
