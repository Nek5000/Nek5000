import os
from subprocess import Popen, PIPE, STDOUT
from pathlib import Path


def build_tools(
    tools_root,
    tools_bin,
    f77=None,
    cc=None,
    bigmem=None,
    targets=("clean", "all"),
    verbose=False,
):

    tools_root = Path(tools_root)

    print("Compiling tools... ")
    print(f'    Using output directory "{tools_bin}"')
    print(f'    Using FC "{f77}"')
    print(f'    Using CC "{cc}"')

    maketools_in = tools_root / "maketools"

    my_env = os.environ.copy()
    if f77:
        my_env["FC"] = f77
    if cc:
        my_env["CC"] = cc
    my_env["bin_nek_tools"] = tools_bin

    if targets[0] == "all":
        targets = os.walk(tools_root).next()[1]
        print(targets)

    for t in targets:
        proc = Popen([maketools_in, t], env=my_env, cwd=tools_root, stderr=STDOUT)
        proc.wait()
        logfile = tools_root / t / "build.log"
        if proc.returncode != 0:
            with open(logfile, "r") as file:
                text = file.read()
            print(text)
            exit(-1)


def build_nek(source_root, usr_file, cwd=None, opts=None, verbose=False):

    if not opts:
        _opts = {}
    else:
        _opts = opts.copy()
    _opts.update(NEK_SOURCE_ROOT=source_root)

    print("Compiling nek5000...")
    print(f'    Using working directory "{cwd}"')
    print(f'    Using .usr file "{usr_file}"')
    for key, val in list(_opts.items()):
        print(f'    Using {key}="{val}"')

    my_env = os.environ.copy()
    if source_root:
        my_env["NEK_SOURCE_ROOT"] = source_root
    if _opts.get("F77"):
        my_env["FC"] = _opts.get("F77")
    if _opts.get("CC"):
        my_env["CC"] = _opts.get("CC")
    if _opts.get("PPLIST"):
        my_env["PPLIST"] = _opts.get("PPLIST")

    makenek_in = Path(source_root) / "bin" / "makenek"
    logfile = Path(cwd) / "build.log"

    proc = Popen([makenek_in, "clean"], cwd=cwd, env=my_env, stdin=PIPE, text=True)
    proc.communicate(input="Y\n")
    proc.wait()

    proc = Popen([makenek_in, usr_file], cwd=cwd, env=my_env, stdin=PIPE, stderr=STDOUT)

    proc.wait()

    if proc.returncode != 0:
        with open(logfile, "r") as file:
            text = file.read()
        print(text)
        exit(-1)
