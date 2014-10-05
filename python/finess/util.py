# coding=utf-8
"""Provides utility functions.
   
   Notably,
   check_all_apps_in_a_dir_compile
"""
from __future__ import absolute_import


def run_command_in_dir(directory, command):
    """Run make_command in directory.
       Returns (command_return_code, stdout_output, stderr_output)
    """
    import os, subprocess
    current_working_directory = os.getcwd()
    os.chdir(directory)
    try:
        p = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        stdout_output, stderr_output = p.communicate()
        command_return_code = p.wait()
    finally:
        os.chdir(current_working_directory)
    return command_return_code, stdout_output, stderr_output

def rebuild_app(app_dir, clean_command = 'make cleanall', build_command = 'make -j'):
    """Rebuild an app.
       Example:
       >>> build_command_return_code, stdout_output_of_build_command, stderr_output_of_build_command = \
           rebuild_app('${FINESS}/apps/2d/euler/test_exact1')
    """
    from os.path import expanduser, expandvars
    app_dir = expanduser(expandvars(app_dir))
    clean_return_code, _, _ = run_command_in_dir(directory = app_dir, command = clean_command.split())
    if clean_return_code != 0:
        raise RuntimeError("Clean failed: " + app_dir)
    return run_command_in_dir(directory = app_dir, command = build_command.split())


def get_app_dirs_in_a_dir(directory):
    """Returns a list of app directories in given directory.

    Criterion: contains both main.cpp and Makefile."""
    import os
    app_dirs = [dir_ for dir_, _, files in os.walk(directory) \
                     if 'main.cpp' in files and \
                        'Makefile' in files]
    return sorted(app_dirs)



def check_all_apps_in_a_dir_compile(directory, show_error_message = True):
    """Rebuild all apps in a directory, and print result to stdout.
       Examples:
            Example 1
            Assume current working directory is $FINESS
            >>> check_all_apps_in_a_dir_compile('apps/2d/advection/')
            Build succeeded: apps/2d/advection/smooth_example
            ========================================
            Build succeeded: apps/2d/advection/rotating
            ========================================
            Build succeeded: apps/2d/advection/sine_example
            ========================================

            Example 2
            Assume current working directory is $FINESS/apps/1d
            >>> check_all_apps_in_a_dir_compile('.')
            Build succeeded: ./advection/c0deriv
            ========================================
            Build succeeded: ./advection/smooth_example
            ========================================
            Build succeeded: ./advection/sq_wave
            ========================================
            Build succeeded: ./buckley_lev/sq_wave
            ========================================
            Build succeeded: ./burgers/shock
            ========================================
            Build succeeded: ./burgers/sine_to_n
            ========================================
            Build succeeded: ./euler/blast
            ========================================
            Build succeeded: ./euler/harten_shock_tube
            ========================================
            Build succeeded: ./euler/shock_entropy
            ========================================
            Build succeeded: ./euler/shock_tube
            ========================================
            Build succeeded: ./euler/test_exact
            ========================================
            Build succeeded: ./mhd/alfven
            ========================================
            Build succeeded: ./mhd/shock_tube
            ========================================
            Build succeeded: ./shallow_water/dam_break
            ========================================

            Example 3
            Assume current working directory is $FINESS.  Also, let's
            say 2d euler apps now fail to compile.
            >>> check_all_apps_in_a_dir_compile("apps/2d/euler")
            Build failed: apps/2d/euler/mach_reflection
            ...
            SetBndValues.cpp: In function ‘void TopFunc(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&)’:
            SetBndValues.cpp:213:27: error: ‘t’ was not declared in this scope
                     if ( x < x0+(20.0*t+1.0)*osq3)
                                       ^
            ...
            make: *** [SetBndValues.o] Error 1
            make: *** Waiting for unfinished jobs....
            ...
            ========================================
            Build failed: apps/2d/euler/radial_shock
            ... more error message ...
            ...
            >>>

    See get_app_dir_in_a_dir for the strategy used to find app
    directories in a given directory.

    Note:
        1. Since 'make -j' was used as the build command, the output
           might be interlaced;
        2. The output message contains only the output of build
           command to stderr, not stdout.
    """
    for app_dir in get_app_dirs_in_a_dir(directory):
        try:
            build_command_return_code, build_command_stdout_output, build_command_stderr_output = \
                rebuild_app(app_dir)
            if build_command_return_code != 0:
                print "Build failed: " + app_dir
                if show_error_message:
                    print build_command_stderr_output
            else:
                print "Build succeeded: " + app_dir
        except RuntimeError as e:
            print e.message

        print "=" * 40

def check_all_apps_compile(show_error_message = True):
    """Rebuild all apps, and print result to stdout."""
    import os
    check_all_apps_in_a_dir_compile(os.environ['FINESS'],
                                    show_error_message)
            
    
    
