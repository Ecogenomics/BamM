#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function

if __name__ == "__main__":

    try:
        import numpy
        print("Numpy version: %s" % numpy.__version__)
        print("Location: %s" % numpy.__file__)
    except ImportError:
        print("Error importing numpy")

    import platform
    import sys
    print(platform.system(), platform.release())
    print((sys.version))
    print((sys.path))
