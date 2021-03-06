function plotfin2(parameters_ini_filename_in)
%PLOTFIN2     Generic 2D plotting routine for FINESS.
%
% PLOTFIN2(outputdir_in) is the main plotting
% routine used in FINESS for all the 2D problems.
%
% In almost every application in FINESS, there is a single script, plotq1.m 
% that this function calls once per frame.  There is a default script located
% here in the library.  That script is where a
% user is expected to modify extra desired plotting stuff, for example setting
% axes, linecolors, exact solutions etc.
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

if( nargin )
  parameters_ini_filename = parameters_ini_filename_in;
else
  parameters_ini_filename = 'parameters.ini';
end

% Read in any user-supplied paramters.
INI = ConvertIniFile2Struct(parameters_ini_filename);

% pull the output directory.
if( isfield( INI.finess, 'output_dir' ) )
  outputdir = INI.finess.output_dir;
else
  outputdir = 'output';
end
INI.finess


% Pull more information from parameters file
ndims = sscanf(INI.finess.ndims, '%d');
if (ndims~=2)
    error(['Incorrect dimension, ndims must be 2. ndims = ',num2str(ndims)]);
end

% Grid and problem information
meqn  = sscanf(INI.finess.meqn, '%d');
maux  = sscanf(INI.finess.maux, '%d');
nplot = sscanf(INI.finess.nout, '%d');
mx    = sscanf(INI.grid.mx, '%d');
my    = sscanf(INI.grid.my, '%d');
xlow  = sscanf(INI.grid.xlow, '%e');
xhigh = sscanf(INI.grid.xhigh, '%e');
ylow  = sscanf(INI.grid.ylow, '%e');
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
      % solution should be found in file
      %     outputdir/q[n1].dat
      fname = [outputdir,'/',num2str(n1+10000),'.dat'];

      % replace the 1000's digit by the letter q
      fname(length(outputdir)+2) = 'q';

      fids = fopen(fname,'r');
      if( fids==-1 )
          error(['File  ',fname,'  not found.']);
      end

      %% Conserved variables -- qsoln
      time   = fscanf(fids,'%e',1);
      q_data = fscanf(fids,'%e',[inf]);
      qsoln = reshape( q_data, mx, my, meqn );
      clear q_data;
      fclose( fids );

      if (maux>0)

        %     outputdir/q[n1].dat
        fname = [outputdir,'/',num2str(n1+10000),'.dat'];

        % replace the 1000's digit by the letter a
        fname(length(outputdir)+2) = 'a';

        fids = fopen(fname,'r');
        if( fids==-1 )
            error(['File  ',fname,'  not found.']);
        end

        %% Aux variables -- aux
        time     = fscanf(fids,'%e',1);
        aux_data = fscanf(fids,'%e',[inf]);
        aux      = reshape( aux_data, mx, my, maux );
        clear aux_data;
        fclose( fids );

      end
      
      % USER SUPPLIED FUNCTION: Plotting function
      plotq2;
    end

disp(' ')

end
