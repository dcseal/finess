function plotfin2_nostop(outputdir_in)
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
% See also: plotq2 (locally defined) and plotfin2

% Read in the parameters file that was used to generate this run
if( ~nargin )
  outputdir_in = 'output/';
end
fids = fopen([[outputdir_in, '/parameters_ini_filename']],'r');
parameters_ini_filename_str = fgetl(fids)
fclose(fids);

% Read in any user-supplied parameters
INI = ConvertIniFile2Struct([[outputdir_in, '/', parameters_ini_filename_str]]);

% Pull the output directory.
%
% TODO - This seems to be a circular argument.  The output_dir is passed in,
% but then it is read in from the parameters file.
if( isfield( INI.finess, 'output_dir' ) )
  outputdir = INI.finess.output_dir;
else
  outputdir = 'output';
end
INI.finess

if( ~strcmp(outputdir, outputdir_in ) )
    disp('Warning: the outputdir read from the parameters file is different than what was passed into this routine.');
    outputdir
    outputdir_in
end


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
m = 1;

for n1=0:nplot

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

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Read in the data from the output folder
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
