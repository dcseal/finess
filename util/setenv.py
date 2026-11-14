# ============================================================================
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import sys
import os
import getopt

def write_environment_variable(csh_handle,bash_handle,var,value):
    csh_handle.write('setenv %s "%s"\n' % (var.upper(),value))
    bash_handle.write('export %s="%s"\n' % (var.upper(),value))

def write_shortcut_alias(csh_handle,bash_handle,var,value):
    csh_handle.write( 'alias %s "%s"\n' % (var,value))
    bash_handle.write('alias %s="%s"\n' % (var,value))


# Paths
finess_path = os.path.abspath(os.curdir)
outfile_base="setenv"
    
# Open output files
csh_file = open(os.path.join(finess_path,".".join((outfile_base,"csh"))),'w')
bash_file = open(os.path.join(finess_path,".".join((outfile_base,"bash"))),'w')

# Write out boiler plate
boiler_plate = ("# FinEss environment settings\n")
csh_file.write(boiler_plate)
bash_file.write(boiler_plate)

# Write out variables
python_path = "${PYTHONPATH}"
matlab_path = "${MATLABPATH}"

finessP = "".join((finess_path,"/python"))
finessM = "".join((finess_path,"/viz/matlab"))

python_path = ":".join((finessP,python_path))
matlab_path = ":".join((finessM,matlab_path))
hdf5_disable_version_check = "1"

print("")
print("The following variables will be set:")
print("  FINESS = %s" % finess_path)
write_environment_variable(csh_file,bash_file,"FINESS",finess_path)

print("  PYTHONPATH = %s" % python_path)
write_environment_variable(csh_file,bash_file,"PYTHONPATH",python_path)

print("  MATLABPATH = %s" % matlab_path)
write_environment_variable(csh_file,bash_file,"MATLABPATH",matlab_path)

print("  HDF5_DISABLE_VERSION_CHECK = %s" % hdf5_disable_version_check)
write_environment_variable(csh_file,bash_file,"HDF5_DISABLE_VERSION_CHECK", hdf5_disable_version_check)

print("")
print("The following aliases will be set:")

# Plotting scripts, all run from an application directory.
# plotfin[n] live in $FINESS/viz/python; plot*_generic use the finess.viz package.
aliases = [
    ("plotfin1",               "$FINESS/viz/python/plotfin1.py"),
    ("plotfin2",               "$FINESS/viz/python/plotfin2.py"),
    ("plot1d_generic",         "$FINESS/viz/plot1d_generic.py"),
    ("plot2d_generic",         "$FINESS/viz/plot2d_generic.py"),
    ("plot2d_generic_contour", "$FINESS/viz/plot2d_generic_contour.py"),
]
for name, script in aliases:
    command = "python3 " + script
    print("  %s = %s" % (name, command))
    write_shortcut_alias(csh_file, bash_file, name, command)
print("")

# Close output files
csh_file.close()
bash_file.close()
