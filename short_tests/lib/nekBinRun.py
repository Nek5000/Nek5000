import os
import sys
from subprocess import call, check_call, PIPE, STDOUT, Popen, CalledProcessError


def run_meshgen(command, stdin, cwd, verbose=False):

    logfile = os.path.join(cwd, "{}.out".format(os.path.basename(command)))
    stdin_bytes = "\n".join(stdin) + "\n"

    print('Running "{}"...'.format(os.path.basename(command)))
    print(f'    Using command "{command}"')
    print(f'    Using stdin "{stdin}"')
    print(f'    Using working directory "{cwd}"')

    try:
        (stdoutdata, stderrdata) = Popen(
            [command], stdin=PIPE, stderr=STDOUT, stdout=PIPE, cwd=cwd, text=True
        ).communicate(stdin_bytes)

        with open(logfile, "w") as file:
            file.writelines(stdoutdata)

        if verbose:
            sys.stdout.write(stdoutdata)

    except (OSError, CalledProcessError) as E:
        # TODO: Change to warnings.warn()
        print(
                'Could not complete {}!  Caught error: "{}".  Check "{}" for details.'.format(
                    command, E, logfile
                )
        )
    else:
        print("Successfully finished {}!".format(os.path.basename(command)))


def run_nek_script(script, rea_file, cwd, log_suffix="", mpi_procs="1"):
    try:
        logs = (
            os.path.join(cwd, "logfile"),
            os.path.join(cwd, f"{rea_file}.log.{mpi_procs}"),
        )

        # Remove old logs
        for l in logs:
            if os.path.exists(l):
                os.remove(l)

        # Running 'script' through shell since it doesn't have a shebang at the top.
        # Need to concatenate args into a string if shell=True
        cmd = " ".join([script, rea_file, str(mpi_procs)])
        print("Running nek5000...")
        print(f'    Using command "{cmd}"')
        print(f'    Using working directory "{cwd}"')
        try:
            # TODO: This doesn't work as intended.  If the nek executable fails, the nek script doesn't return the error.
            # Check doxygen to see what exit values there are (some succesful exit values there are!)
            check_call(cmd, cwd=cwd, shell=True)
        except Exception as E:
            # TODO: Change to warnings.warn()
            print(f"Could not successfully run nek5000! Caught error: {E}")
        else:
            # print('Successfully ran nek5000!')
            print("Finished running nek5000!")

        # Rename logs
        if log_suffix:
            for l in logs:
                os.rename(l, l + log_suffix)

    # This are expected exceptions if 'check_call' or 'os.rename' fail.
    # We issue a warning, not error, so subsequent tests can continue
    except (OSError, CalledProcessError) as E:
        # TODO: Change to warnings.warn()
        print(
                'Could not complete command: "{}": {}'.format(
                    " ".join([script, rea_file, mpi_procs]), E
                )
        )


def run_nek(
    cwd, rea_file, ifmpi, log_suffix="", n_procs=1, step_limit=None, verbose=False
):
    # Paths to executables, files
    nek5000 = os.path.join(cwd, "nek5000")
    logfile = os.path.join(cwd, f"{rea_file}.log.{n_procs}{log_suffix}")
    session_name = os.path.join(cwd, "SESSION.NAME")
    ioinfo = os.path.join(cwd, "ioinfo")
    if ifmpi:
        command = ["mpiexec", "-np", str(n_procs), nek5000]
    else:
        command = [nek5000]

    print("Running nek5000...")
    print('    Using command "{}"'.format(" ".join(command)))
    print(f'    Using working directory "{cwd}"')
    print(f'    Using .rea file "{rea_file}"')

    # Any error here is unexepected
    try:
        with open(session_name, "w") as f:
            f.writelines(
                [
                    "{}\n".format(1),
                    f"{rea_file}\n",
                    "{}\n".format(cwd + "/"),
                ]
            )

        if step_limit:
            with open(ioinfo, "w") as f:
                f.writelines([f"-{step_limit}"])

        if verbose:
            with open(logfile, "w") as f:
                proc = Popen(command, cwd=cwd, stderr=STDOUT, stdout=PIPE, text=True)
                for line in proc.stdout:
                    sys.stdout.write(line)
                    f.write(line)
        else:
            with open(logfile, "w") as f:
                call(command, cwd=cwd, stdout=f)

    except Exception as E:
        # TODO: Change to warnings.warn()
        print(f"Could not successfully run nek5000! Caught error: {E}")
    else:
        print("Finished running nek5000!")


def run_neknek(
    cwd,
    inside,
    np_inside,
    outside,
    np_outside,
    coupled=True,
    step_limit=None,
    log_suffix="",
    verbose=False,
):

    # Paths to executables, files
    nek5000 = os.path.join(cwd, "nek5000")
    logfile = os.path.join(
        cwd,
        "{inside}{np_in}.{outside}{np_out}.log{sfx}".format(
            inside=inside,
            outside=outside,
            np_in=np_inside,
            np_out=np_outside,
            sfx=log_suffix,
        ),
    )

    ifcoupled = "F"
    if coupled:
        ifcoupled = "T"

    inside_log = os.path.join(cwd, f"{inside}.log")
    inside_his = os.path.join(cwd, f"{inside}.his")

    outside_log = os.path.join(cwd, f"{outside}.log")
    outside_his = os.path.join(cwd, f"{outside}.his")

    session_name = os.path.join(cwd, "SESSION.NAME")
    ioinfo = os.path.join(cwd, "ioinfo")

    command = ["mpiexec", "-np", str(int(np_inside) + int(np_outside)), nek5000]

    print("Running nek5000...")
    print('    Using command "{}"'.format(" ".join(command)))
    print(f'    Using working directory "{cwd}"')
    print(f'    Using .rea files "{inside}", "{outside}"')

    # Any error here is unexpected
    try:

        # Create SESSION.NAME
        with open(session_name, "w") as f:
            f.writelines(
                [
                    "{}\n".format(2),
                    f"{ifcoupled}\n",
                    f"{inside}\n",
                    f"{cwd}\n",
                    f"{np_inside}\n",
                    f"{outside}\n",
                    f"{cwd}\n",
                    f"{np_outside}\n",
                ]
            )

        # Write step limit
        if step_limit:
            with open(ioinfo, "w") as file:
                file.writelines([f"-{step_limit}"])

        if verbose:
            with open(logfile, "w") as f:
                proc = Popen(command, cwd=cwd, stderr=STDOUT, stdout=PIPE, text=True)
                for line in proc.stdout:
                    sys.stdout.write(line)
                    f.write(line)
        else:
            with open(logfile, "w") as f:
                call(command, cwd=cwd, stdout=f)

    except Exception as E:
        # TODO: Change to warnings.warn()
        print(f"Could not successfully run nek5000! Caught error: {E}")
    else:
        print("Finished running nek5000!")


def mvn(src_prefix, dst_prefix, cwd):
    exts = (".box", ".rea", ".usr", ".map", ".sep", ".re2")
    print("Running mvn...")
    print(f'    Using working directory "{cwd}"')
    for x in exts:
        src = os.path.join(cwd, src_prefix + x)
        dst = os.path.join(cwd, dst_prefix + x)
        try:
            os.rename(src, dst)
        except OSError as E:
            # TODO: Change to warnings.warn()
            print(f"    Could not move {src} to {dst}: {E}")
        else:
            print(f"    Successfully moved {src} to {dst}")
    print("Finished running mvn!")
