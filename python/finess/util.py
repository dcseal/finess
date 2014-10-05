"""Provides utility functions.
"""
from __future__ import absolute_import


def run_command_in_dir(directory, command):
    """Run make_command in directory.
       Returns (command_return_code, stdin_output, stderr_output)
    """
    import os, subprocess
    current_working_directory = os.getcwd()
    os.chdir(directory)
    try:
        p = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        stdin_output, stderr_output = p.communicate()
        command_return_code = p.wait()
    finally:
        os.chdir(current_working_directory)
    return command_return_code, stdin_output, stderr_output

def rebuild_app(app_dir, clean_command = 'make cleanall', build_command = 'make -j'):
    """Rebuild an app.
       Example:
       >>> build_command_return_code, stdin_output_of_build_command, stderr_output_of_build_command = \
           rebuild_app('${FINESS}/apps/2d/euler/test_exact1')
    """
    from os.path import expanduser, expandvars
    app_dir = expanduser(expandvars(app_dir))
    clean_return_code, _, _ = run_command_in_dir(directory = app_dir, command = clean_command.split())
    if clean_return_code != 0:
        raise RuntimeError("Clean failed: " + app_dir)
    return run_command_in_dir(directory = app_dir, command = build_command.split())


def get_app_dirs_in_a_dir(directory):
    "Returns a list of strings.  Criterion: contains both generate_iniparams.py and Makefile."
    import os
    app_dirs = [dir_ for dir_, _, files in os.walk(directory) \
                     if 'generate_iniparams.py' in files and \
                        'Makefile' in files]
    return app_dirs



def check_all_apps_in_a_dir_compile(directory, show_error_message = True):
    """Rebuild all apps in a directory, and print result to stdin."""
    for app_dir in get_app_dirs_in_a_dir(directory):
        try:
            build_command_return_code, build_command_stdin_output, build_command_stderr_output = \
                rebuild_app(app_dir)
            if build_command_return_code != 0:
                print "Build failed: " + app_dir
                if show_error_message:
                    print build_command_stderr_output
            else:
                print "Build succeeded: " + app_dir
        except RuntimeError as e:
            print e.value

        print "=" * 40

def check_all_apps_compile(show_error_message = True):
    """Rebuild all apps, and print result to stdin."""
    import os
    check_all_apps_in_a_dir_compile(os.environ['FINESS'],
                                    show_error_message)
            
    
    
