function plotfin1(outputdir_in)
%PLOTFIN1    Generic 1D plotting routine for FINESS.
%
% PLOTFIN1( outputdir_in ) is the main MATLAB plotting routine used in FINESS
% for all the 1D problems.
%
% outputdir_in = location for output directory.  default = 'output'
%
% In almost every application in FINESS, there is a single script, plotq1.m 
% that this function calls once per frame.  There is a default script located
% here in the library.  That script is where a
% user is expected to modify extra desired plotting stuff, for example setting
% axes, linecolors, exact solutions etc.
%
% This script was borrowed and modified from DoGPack (dogpack-code.org).
%
% Input parameters:
%
%   outputdir_in - string identifying the location of the
%   output directory.  The parameters file will then be read in by examining
%   the parameters_ini_filename that FINESS creates when it runs the code.
%   Default = 'output'.
%
% Output parameters:
%
%    None.
%
% See also: plotq1, plotfin2, for the other dimensional plotters.

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
if( isfield( INI.finess, 'output_dir' ) )
    outputdir = INI.finess.output_dir;
    if( ~strcmp(outputdir, outputdir_in ) )
        disp('Warning: the outputdir read from the parameters file is different than what was passed into this routine.');
        disp('    Using the file that is passed into this function.');
        outputdir
        outputdir_in
    end
    outputdir = outputdir_in;
    INI.finess.output_dir = outputdir;
else
  outputdir = 'output';
end
INI.finess

% Pull more information from parameters file
ndims = sscanf(INI.finess.ndims, '%d');
if (ndims~=1)
    error(['Incorrect dimension, ndims must be 1. ndims = ',num2str(ndims)]);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%  Parse QHELP.DAT (or parameters.ini file)
% %%  THIS HAS BEEN DECPRECATED IN FAVOR OF USING DATA FROM PARAMETERS FILE.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fids  = fopen([outputdir,'/qhelp.dat'],'r');
% if( fids==-1 )
%     error(['File  ',outputdir,'/qhelp.dat  not found.']);
% end
% ndims = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
% if( ndims~=1 )
%     error(['Incorrect dimension, ndims must be 1. ndims = ',num2str(ndims)]);
% end
% meqn    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
% maux    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
% nplot   = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
% %
% mx      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
% xlow    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
% xhigh   = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
% fclose(fids);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Grid and problem information
  meqn  = sscanf(INI.finess.meqn, '%d');
  maux  = sscanf(INI.finess.maux, '%d');
  nplot = sscanf(INI.finess.nout, '%d');
  mx    = sscanf(INI.grid.mx, '%d');
  xlow  = sscanf(INI.grid.xlow, '%e');
  xhigh = sscanf(INI.grid.xhigh, '%e');

  % (Uniform) Grid information
  dx = (xhigh-xlow)/mx;
  xc = transpose(linspace(xlow+dx/2,xhigh-dx/2,mx));

  % q - flag for determining whether or not to quit
  q=-1;
  m = input([ 'Which component of q do you want to plot ( 1 - ',num2str(meqn),...
              ' ) ? ']);
    disp(' ')
  if isempty(m)
    m=1;
  end

  nf = 0;       % Frame number
  n1 = -1;      % Frame number

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
      if fids==-1
          error(['File  ',fname,'  not found.']);
      end
 

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Read in the data from the output folder
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %% Conserved variables -- qsoln
      time  = fscanf(fids,'%e', 1);
      qtmp  = fscanf(fids,'%e', [1,inf]);
      qsoln = reshape( qtmp, mx, meqn );
      clear qtmp;

      %% Aux variables -- aux
      if( maux > 0 )
          fname(length(outputdir)+2) = 'a';
          fids = fopen(fname,'r');
          if( fids == -1 )
              error(['File  ',fname,'  not found.']);
          end
          time = fscanf(fids,'%e',1);
          atmp = fscanf(fids,'%e',[1,inf]);
          aux  = reshape( atmp, mx, maux );
          clear atmp;

      end
  
      % USER SUPPLIED FUNCTION
      plotq1;        
      end
  
   end
   disp(' ')

end
