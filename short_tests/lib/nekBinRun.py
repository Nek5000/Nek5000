import os
import sys
from warnings import warn
from subprocess import call, check_call, PIPE, STDOUT, Popen, CalledProcessError
from pathlib import Path


def run_meshgen(command, stdin, cwd, verbose=False):

    base_command = Path(command).name
    logfile = Path(cwd) / f"{base_command}.out"

    print(
        f'Running "{base_command}"...\n'
        f'    Using command "{command}"\n'
        f'    Using stdin "{stdin}"\n'
        f'    Using working directory "{cwd}"'
    )
    stdin = "\n".join(stdin) + "\n"

    try:
        (stdoutdata, stderrdata) = Popen(
            [command], stdin=PIPE, stderr=STDOUT, stdout=PIPE, cwd=cwd, text=True
        ).communicate(stdin)

        with open(logfile, "w") as file:
            file.writelines(stdoutdata)

        if verbose:
            sys.stdout.write(stdoutdata)

    except (OSError, CalledProcessError) as error:
        warn(
            f'Could not complete {command}!  Caught error: "{error}".  '
            f'Check "{logfile}" for details.'
        )
        raise
    else:
        print(f"Successfully finished {base_command}!")


def run_nek_script(script, rea_file, cwd, log_suffix="", mpi_procs="1"):
    cwd = Path(cwd)
    try:
        logs = (
            cwd / "logfile",
            cwd / f"{rea_file}.log.{mpi_procs}",
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
        except Exception as error:
            warn(f"Could not successfully run nek5000! Caught error: {error}")
        else:
            # print('Successfully ran nek5000!')
            print("Finished running nek5000!")

        # Rename logs
        if log_suffix:
            for l in logs:
                os.rename(l, l + log_suffix)

    # This are expected exceptions if 'check_call' or 'os.rename' fail.
    # We issue a warning, not error, so subsequent tests can continue
    except (OSError, CalledProcessError) as error:
        warn(
            'Could not complete command: "{}": {}'.format(
                " ".join([script, rea_file, mpi_procs]), error
            )
        )


def run_nek(
    cwd, rea_file, ifmpi, log_suffix="", n_procs=1, step_limit=None, verbose=False
):
    # Paths to executables, files
    cwd = Path(cwd)
    nek5000 = str(cwd / "nek5000")
    logfile = cwd / f"{rea_file}.log.{n_procs}{log_suffix}"
    session_name = cwd / "SESSION.NAME"
    ioinfo = cwd / "ioinfo"
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
        with open(session_name, "w") as file:
            file.writelines(
                [
                    "1\n",
                    f"{rea_file}\n",
                    f"{cwd}/\n",
                ]
            )

        if step_limit:
            with open(ioinfo, "w") as file:
                file.writelines([f"-{step_limit}"])

        if verbose:
            with open(logfile, "w") as file:
                proc = Popen(command, cwd=cwd, stderr=STDOUT, stdout=PIPE, text=True)
                for line in proc.stdout:
                    sys.stdout.write(line)
                    file.write(line)
        else:
            with open(logfile, "w") as file:
                call(command, cwd=cwd, stdout=file)

    except Exception as error:
        warn(f"Could not successfully run nek5000! Caught error: {error}")
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
    cwd = Path(cwd)
    nek5000 = str(cwd / "nek5000")
    logfile = cwd / f"{inside}{np_inside}.{outside}{np_outside}.log{log_suffix}"

    ifcoupled = "F"
    if coupled:
        ifcoupled = "T"

    session_name = cwd / "SESSION.NAME"
    ioinfo = cwd / "ioinfo"

    command = ["mpiexec", "-np", str(int(np_inside) + int(np_outside)), nek5000]

    print("Running nek5000...")
    print('    Using command "{}"'.format(" ".join(command)))
    print(f'    Using working directory "{cwd}"')
    print(f'    Using .rea files "{inside}", "{outside}"')

    # Any error here is unexpected
    try:

        # Create SESSION.NAME
        with open(session_name, "w") as file:
            file.writelines(
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
            with open(logfile, "w") as file:
                proc = Popen(command, cwd=cwd, stderr=STDOUT, stdout=PIPE, text=True)
                for line in proc.stdout:
                    sys.stdout.write(line)
                    file.write(line)
        else:
            with open(logfile, "w") as file:
                call(command, cwd=cwd, stdout=file)

    except Exception as error:
        warn(f"Could not successfully run nek5000! Caught error: {error}")
    else:
        print("Finished running nek5000!")

def run_neknekn(cwd, in1, np_in1, in2, np_in2, outside, np_outside, coupled=True, step_limit=None, log_suffix='', verbose=False):

    # Paths to executables, files
    nek5000 = os.path.join(cwd, 'nek5000')
    logfile      = os.path.join(cwd, '{in1}{np_in1}.{in2}{np_in2}.{outside}{np_out}.log{sfx}'.format(
        in1 = in1,
        in2 = in2,
        outside = outside,
        np_in1 = np_in1,
        np_in2 = np_in2,
        np_out = np_outside,
        sfx = log_suffix
    ))

   
    ifcoupled = 'F'
    if coupled :
        ifcoupled = 'T'

    in1_log = os.path.join(cwd, '{0}.log'.format(in1))
    in1_his = os.path.join(cwd, '{0}.his'.format(in1))

    in2_log = os.path.join(cwd, '{0}.log'.format(in2))
    in2_his = os.path.join(cwd, '{0}.his'.format(in2))

    outside_log = os.path.join(cwd, '{0}.log'.format(outside))
    outside_his = os.path.join(cwd, '{0}.his'.format(outside))

    session_name = os.path.join(cwd, 'SESSION.NAME')
    ioinfo       = os.path.join(cwd, 'ioinfo')

    command = ['mpiexec', '-np', str(int(np_in1) + int(np_in2) + int(np_outside)), nek5000]

    print("Running nek5000...")
    print('    Using command "{0}"'.format(' '.join(command)))
    print('    Using working directory "{0}"'.format(cwd))
    print('    Using .rea files "{0}", "{1}"'.format(in1, in2, outside))

    # Any error here is unexpected
    try:

        # Create SESSION.NAME
        with open(session_name, 'w') as f:
            f.writelines([
                "{0}\n".format(3),
                "{0}\n".format(ifcoupled),
                "{0}\n".format(in1),
                "{0}\n".format(cwd),
                "{0}\n".format(np_in1),
                "{0}\n".format(in2),
                "{0}\n".format(cwd),
                "{0}\n".format(np_in2),
                "{0}\n".format(outside),
                "{0}\n".format(cwd),
                "{0}\n".format(np_outside),
            ])

        # Write step limit
        if step_limit:
            with open(ioinfo, 'w') as f:
                f.writelines(['-{0}'.format(step_limit)])

        if verbose:
            with open(logfile, 'w') as f:
                proc =Popen(command, cwd=cwd, stderr=STDOUT, stdout=PIPE)
                for line in proc.stdout:
                    sys.stdout.write(line)
                    f.write(line)
        else:
            with open(logfile, 'w') as f:
                call(command, cwd=cwd, stdout=f)

    except Exception as E:
        # TODO: Change to warnings.warn()
        print('Could not successfully run nek5000! Caught error: {0}'.format(E))
    else:
        print('Finished running nek5000!')

def mvn(src_prefix, dst_prefix, cwd):
    exts = (".box", ".rea", ".usr", ".map", ".sep", ".re2")
    print("Running mvn...")
    print(f'    Using working directory "{cwd}"')
    cwd = Path(cwd)
    for x in exts:
        src = cwd / (src_prefix + x)
        dst = cwd / (dst_prefix + x)
        try:
            src.rename(dst)
        except OSError as error:
            warn(f"    Could not move {src} to {dst}: {error}")
        else:
            print(f"    Successfully moved {src} to {dst}")
    print("Finished running mvn!")
