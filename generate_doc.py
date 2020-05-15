#!/usr/bin/env python
"""
This program generates the html documentation
"""

import subprocess
import sys
import os


def generate_documentation() -> bool:
    """
    Generate html documentation.
    Returns true if the doc has been correctly generated. False otherwise
    """
    # Go to doc repository and launch the command to generate doc
    os.chdir(f"xfv/doc")
    print(f"Working dir is " + os.getcwd())
    launch_cmd = f"make html"
    print(f"Executing command: {launch_cmd}")
    try:
        subprocess.check_output(launch_cmd.split())
    except subprocess.CalledProcessError:
        print(f"The command {launch_cmd} has failed!", file=sys.stderr)
        return False
    return True


if __name__ == "__main__":
    if not generate_documentation():
        sys.exit(1)
    sys.exit(0)
