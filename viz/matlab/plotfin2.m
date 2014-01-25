function plotfin2(outputdir_in)
%PLOTFIN2     Generic 2D plotting routine for FINESS.
%
% PLOTFIN2(outputdir_in) is the main plotting
% routine used in FINESS for all the 2D problems.
%
% Input parameters:
%
%   outputdir_in - string identifying the output directory.  This can be
%   relative or an absolute pathname.  Default = 'output'.
%
% Output parameters:
%
%    None.
%
% See also: plotfin1


  % Parse the input parameters (set default values for points_per_dir,
  % point_type and output folder name:
  global outputdir
  if(nargin<1)
    outputdir = 'output';
  else
    outputdir = outputdir_in;
  end

  % hard-coded values here.  TODO - these should be removed.
  points_per_dir = 1;
  point_type     = 1;
  kmax           = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FIND OUT IF CARTESIAN OR UNSTRUCTURED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fids  = fopen([outputdir,'/qhelp.dat'],'r');
if fids==-1
  error(['File  ',outputdir,'/qhelp.dat  not found.']);
end
ndims = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
if (ndims~=2)
    error(['Incorrect dimension, ndims must be 2. ndims = ',num2str(ndims)]);
end
GridType = fscanf(fids,'%s',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
GridTmp = GridType(1:9);
if (GridTmp=='Cartesian')
  GridType='Cartesian   ';
elseif (GridTmp=='Unstructu')
  GridType='Unstructured';
else
  error(['Incorrect GridType, must be Cartesian or Unstructured. GridType = ',GridType]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');
disp(['        GridType = ',GridType]);
disp(['  points_per_dir = ',num2str(points_per_dir)]);
disp(['      point_type = ',num2str(point_type)]);
disp(['       outputdir = ',outputdir]);
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CARTESIAN plotting routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (GridType=='Cartesian   ')
    
  meqn    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  maux    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  nplot   = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  meth1   = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  datafmt = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  %
  mx      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  my      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  xlow    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  xhigh   = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  ylow    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  yhigh   = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  fclose(fids);
  
  % (Uniform) Grid information
  dx = (xhigh-xlow)/mx;
  dy = (yhigh-ylow)/my;
  xc = linspace(xlow+dx/2,xhigh-dx/2,mx);
  yc = linspace(ylow+dy/2,yhigh-dy/2,my);
  [xc,yc]=meshgrid(xc,yc);
  xc = xc';
  yc = yc';    

  % Default plotting limits (if none are specified in the application)
  xeps = max(0.015*(xhigh-xlow),0.015*(yhigh-ylow));
  yeps = xeps;

  % flag to determine whether or not to quit the program
  q=-1;

  m = input([ 'Which component of q do you want to plot ( 1 - ',num2str(meqn),...
              ' ) ? ']);
  disp(' ')
  if isempty(m)
    m=1;
  end

  n  = 0;
  nf = 0;
  n1 = -1;

  while(nf~=-1)
    nf  = input([ ' Plot which frame ( 0 - ',num2str(nplot),...
                  ' ) [type -1 or q to quit] ? ']);
    if isempty(nf)
      n1 = n1 + 1;
      nf = 0;
    else
      n1 = nf;
    end
    if n1> nplot
      disp(' ');
      disp(' End of plots ');
      disp(' ');
      n1 = nplot;
    end
    if (nf~=-1)
      
      %% Solution -- q
      [q_data, time] = read_state2_cart(datafmt, outputdir, n1, 'q', ...
                                       mx, my, meqn, kmax, ...
                                       1:meqn);
      qsoln = reshape( q_data, mx, my, meqn );
      clear q_data;

      qaug = zeros(mx+1,my+1,meqn);
      qaug(1:mx,1:my,1:meqn) = qsoln;

      if (maux>0)
        %% Aux variables -- aux
        [a_data,time] = read_state2_cart(datafmt, outputdir, n1, 'a', ...
                                         mx, my, maux, kmax, 1:maux);
        aux = reshape( a_data, mx, my, maux );
        clear a_data;
        aux_aug = zeros(mx+1,my+1,maux);
        aux_aug(1:mx,1:my,1:maux) = aux;
      end
      
      % USER SUPPLIED FUNCTION: Plotting function
      plotq2_cart;
    end

  end
  disp(' ')

else
  disp(' ');
  disp([' Error in plotdog2.m: GridType = ',GridType,' is not ' ...
                      'supported.']);
  disp(' ');
end

end
