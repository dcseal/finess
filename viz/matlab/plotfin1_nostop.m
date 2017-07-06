function plotfin1_nostop( parameters_ini_filename_in )
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
%   parameters_ini_filename_in - string identifying the location of the
%   parameters file used to construct this solution.  This can be
%   relative or an absolute pathname.  Default = 'parameters.ini'.
%
% Output parameters:
%
%    None.
%
% See also: plotq1, plotfin2, for the other dimensional plotters.

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

  % Pull more information from parameters file
  ndims = sscanf(INI.finess.ndims, '%d');
  if (ndims~=1)
      error(['Incorrect dimension, ndims must be 1. ndims = ',num2str(ndims)]);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  Parse QHELP.DAT (or parameters.ini file)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fids  = fopen([outputdir,'/qhelp.dat'],'r');
  if( fids==-1 )
      error(['File  ',outputdir,'/qhelp.dat  not found.']);
  end
  ndims = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  if( ndims~=1 )
      error(['Incorrect dimension, ndims must be 1. ndims = ',num2str(ndims)]);
  end
  meqn    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  maux    = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  nplot   = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  %
  mx      = fscanf(fids,'%d',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  xlow    = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  xhigh   = fscanf(fids,'%e',1); fscanf(fids,'%s',1); fscanf(fids,'%s',1);
  fclose(fids);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Grid information
  dx = (xhigh-xlow)/mx;                                 % cell spacing
  xc = transpose(linspace(xlow+dx/2, xhigh-dx/2, mx));  % cell centers

  m = 1;

  for n1=0:nplot

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
      time  = fscanf(fids,'%e', 1)
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
    disp(' ')

end
