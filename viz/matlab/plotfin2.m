function plotfin2(parameters_ini_filename_in)
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
  if(nargin<1)
    parameters_ini_filename = 'parameters.ini';
  else
    parameters_ini_filename = parameters_ini_filename_in;
  end
 
  global INI;
  INI = ini2struct(parameters_ini_filename);
  
  global outputdir;
  if(isfield(INI.finess, 'output_dir'))
    outputdir = INI.finess.output_dir;
  else
    outputdir = 'output';
  end

  % hard-coded values here.  TODO - these should be removed.
  points_per_dir = 1;
  point_type     = 1;
  kmax           = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FIND OUT IF CARTESIAN OR UNSTRUCTURED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ndims = sscanf(INI.finess.ndims, '%d');

if (ndims~=2)
    error(['Incorrect dimension, ndims must be 2. ndims = ',num2str(ndims)]);
end
GridType = 'Cartesian'
GridTmp = GridType(1:9);
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
if (GridType=='Cartesian')

  meqn = sscanf(INI.finess.meqn, '%d');
  maux = sscanf(INI.finess.maux, '%d');
  nplot = sscanf(INI.finess.nout, '%d');
%  meth1 = seems unused
  datafmt = 1
  mx = sscanf(INI.grid.mx, '%d');
  my = sscanf(INI.grid.my, '%d');
  xlow = sscanf(INI.grid.xlow, '%e');
  xhigh = sscanf(INI.grid.xhigh, '%e');
  ylow = sscanf(INI.grid.ylow, '%e');
  yhigh = sscanf(INI.grid.yhigh, '%e');


 
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
