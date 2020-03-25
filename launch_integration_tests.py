#!/usr/bin/env python
"""
This program launches the integration tests suite
"""
from contextlib import contextmanager
import os
from pathlib import Path
import subprocess
import sys
from typing import List


@contextmanager
def down_to_directory(directory: Path) -> None:
    """
    A context manager that descends into the directory and returns to the original one
    """
    current_path = os.getcwd()
    os.chdir(directory.as_posix())
    yield
    os.chdir(current_path)


def get_integration_tests() -> List[Path]:
    """
    Returns the list of paths toward integration tests directories
    """
    integration_tests = []
    integration_root = Path("./xfv/tests/integration")
    for root, _, files in os.walk(integration_root.as_posix()):
        if 'XDATA.xml' in files:
            integration_tests.append(Path(root))
    return integration_tests


def launch_integration_test() -> bool:
    """
    Launch the integration test and return True if everything is ok. False otherwise.
    """
    launch_cmd = "XtendedFiniteVolume.py $(pwd)"
    print(f"Executing command: {launch_cmd}")
    try:
        subprocess.check_output(launch_cmd.split())
    except subprocess.CalledProcessError:
        print(f"The command {launch_cmd} has failed!", file=sys.stderr)
        return False
    compar_cmd = "h5diff reference.hdf5 all_fields.hdf5"
    print(f"Executing command: {compar_cmd}")
    try:
        subprocess.check_output(compar_cmd.split())
    except subprocess.CalledProcessError:
        print("Differences are detected between reference and obtained results!")
        return False
    return True


if __name__ == "__main__":
    for integration_test in get_integration_tests():
        with down_to_directory(integration_test):
            print(f"Launching the integration test {integration_test.as_posix()}")
            if not launch_integration_test():
                sys.exit(1)
    sys.exit(0)
