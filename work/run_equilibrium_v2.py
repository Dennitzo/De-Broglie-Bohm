"""Bootstrap: load double_slit_equilibrium_v2.py as a module and call main().

The isolated research runner adds neither cwd nor the script directory to
sys.path. This wrapper only adjusts sys.path and delegates to the unmodified
entrypoint.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import importlib.util  # noqa: E402

spec = importlib.util.spec_from_file_location(
    "double_slit_equilibrium_v2",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "double_slit_equilibrium_v2.py"),
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.main()
