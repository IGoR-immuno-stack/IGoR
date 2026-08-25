python
import sys
import os
import glob
# Adapted from: https://sourceware.org/gdb/wiki/STLSupport

# Find pixi environment relative to cwd
cwd = os.getcwd()
patterns = [
    '.pixi/envs/default/share/gcc-*/python',
    '.pixi/envs/*/share/gcc-*/python',
]

for pattern in patterns:
    matches = glob.glob(os.path.join(cwd, pattern))
    if matches:
        for path in matches:
            if os.path.exists(os.path.join(path, 'libstdcxx/v6/printers.py')):
                sys.path.insert(0, path)
                from libstdcxx.v6.printers import register_libstdcxx_printers
                register_libstdcxx_printers(None)
                print(f"SUCCESS: Registered STL printers from {path}")
                break
        break
else:
    print("WARNING: Could not find libstdcxx printers. STL pretty-printing disabled.")
end