#! /usr/bin/env python
"""Git pre-commit hook to run pylint on python files.

Slight variation of: 
    https://gist.github.com/nivbend/7e0e306a98138916b3c9#file-run_pylint-py
    
Copy it under .git/hooks/pre-commit at the root of your project.
"""

from __future__ import print_function
from subprocess import check_output, CalledProcessError
from sys import stderr
from os.path import isfile
from os import getenv

(SUCCESS,
 GIT_DIFF_ERROR,
 PYLINT_ERRORS) = range(3)

PYLINTRC = ".pylintrc"

def _print_error(message):
    """Print an error message to stderr."""
    print(message, file = stderr)

def _is_python_script(filename):
    """Return true for *.py files and python scripts ("#! /path/to/python")."""
    if not isfile(filename):
        return False

    if not filename.endswith(".py"):
        try:
            with open(filename, "rb") as contents:
                first_line = contents.readline()
        except OSError:
            return False

        # Check shebang.
        if not (first_line.startswith("#!") and "python" in first_line):
            return False

    return True


def get_git_cmd():
    """Return the git command to detect the changes between branch and master"""
    # is_travis_env = getenv("$TRAVIS_BRANCH") is not None
    # if is_travis_env:
    cmd = "git diff --name-only --diff-filter=AM master-new"
    # else:
    #     cmd = "git diff --staged --name-only HEAD"
    print ("cmd is {:s}".format(cmd))
    return cmd.split()


def run():
    """Verify changed python files using pylint."""
    # Get all changed files' paths.
    try:
        changed_files = check_output(get_git_cmd())
    except CalledProcessError:
        _print_error("Couldn't get list of changed files")
        return GIT_DIFF_ERROR

    # Limit checks to python scripts only.
    changed_files = [
        filename for filename
        in changed_files.splitlines()
        if _is_python_script(filename)]

    if changed_files:
        cmd = ["pylint",]
        if isfile(PYLINTRC):
            cmd += ["--rcfile={}".format(PYLINTRC)]
        try:
            check_output(cmd + changed_files)
        except CalledProcessError as error:
            _print_error(error.output)
            _print_error("pylint returned errors, aborting commit.")
            return PYLINT_ERRORS

    return SUCCESS

if __name__ == "__main__":
    exit(run())
