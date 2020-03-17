# -*- coding: iso-8859-1 -*-
"""
Different utilities for testing
"""
from contextlib import contextmanager
from StringIO import StringIO
import sys

@contextmanager
def captured_output():
    """
    A context manager for capturing output (by Jonathon Reinhart)
    """
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
         sys.stdout, sys.stderr = old_out, old_err

