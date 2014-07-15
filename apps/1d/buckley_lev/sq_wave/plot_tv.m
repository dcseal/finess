function tv = plot_tv(outputdir)

    format long e

    % INFO --------------------------------------------------------------
    disp(' ')
    disp('TOTAL VARIATION PLOTTER')
    disp(' ')

    % meqn and nplot
    fids  = fopen(strcat(outputdir, '/qhelp.dat'),'r');
    if fids==-1
        error(['File  ',outputdir,'/qhelp.dat  not found.']);
    end

    nplot = fscanf(fids,'%d',1);
    meqn  = fscanf(fids,'%d',1);
    fclose(fids);

    m = input([ 'Which component of q do you want to plot ( 1 - ',num2str(meqn),...
        ' ) ? ']);
    disp(' ')
    if isempty(m)
      m=1;
    end

    % INFO --------------------------------------------------------------

    fid = fopen( strcat(outputdir, '/total_variation.dat' ) );
    consv = fscanf(fid,'%e',[meqn+1 inf]);
    status = fclose(fid);

    t = consv(1,:)';
    qc = consv(m+1,:)';
    nt = length(t);

    for i=1:nt
      dqc(i) = qc(i)-qc(1);
    end

    figure(2)
    clf
    plot(t,dqc,'b-');
    set(gca,'fontsize',16);
    t1=title([ 'Total variation of q(',num2str(m),') vs. time']);
    set(t1,'fontsize',16);
    axis on; box on; grid off;
    hold off

end
