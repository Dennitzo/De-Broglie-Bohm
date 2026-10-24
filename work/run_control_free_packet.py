"""Bootstrap: load control_free_packet.py as a module and call main().

The isolated research runner adds neither cwd nor the script directory to
sys.path. This wrapper only adjusts sys.path and delegates to the unmodified
entrypoint; progress prints help locate failures.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import importlib.util  # noqa: E402

spec = importlib.util.spec_from_file_location(
    "control_free_packet", os.path.join(os.path.dirname(os.path.abspath(__file__)), "control_free_packet.py")
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
print("MODULE LOADED")
module.main()
print("MAIN DONE")
