"""CLI for logs."""
import sys
import textwrap

from .. import utils as cli_utils
from ... import tasks


def stream_task_logs_tail(args):
    task_id = args.mode.lower()
    task = tasks.Task(task_id)

    filename = "stdout.txt" if args.stdout else "stderr.txt"
    files = ["stdout.txt", "stderr.txt"
            ] if args.stdout == args.stderr else [filename]

    return task.tail_files(files, 10, True, sys.stdout)


def register(parser):
    cli_utils.show_help_msg(parser)

    parser.add_argument(
        "mode",
        type=str,
        nargs="?",
        default="SUBMITTED",
        help=
        ("Mode of log retrieval.\n"
         "Use 'SUBMITTED' for the last submitted task, 'SUBMITTED-1' for the \n"
         "second last submitted task, and so on. 'STARTED' and 'STARTED-n'\n"
         "follow the same pattern for started tasks. Or, use a specific\n"
         "task ID to retrieve logs for a particular task."))

    parser.add_argument("--stdout",
                        action="store_true",
                        help="Stream the task's standard output (STDOUT) in " \
                             "real time.")
    parser.add_argument("--stderr",
                        action="store_true",
                        help="Stream the task's standard error (STDERR) in "
                        "real time.")
    parser.add_argument("--no-color",
                        action="store_true",
                        help="Disable colored output.")
    parser.add_argument(
        "--wait",
        "-w",
        action="store_true",
        help=("Wait for the task to start running before streaming logs.\n"
              "Without this flag, logs are streamed immediately, or an error\n"
              "is returned if the task is not running."))
    # Register function to call when this subcommand is used

    parser.epilog = textwrap.dedent("""\
        examples:

            # Real-time streaming of standard output from an OpenFOAM simulation
            $ inductiva logs f0bnqgf4fcr4asgi4e21tcsqa
            Websocket connected
            Opening socket connection to logs of task f0bnqgf4fcr4asgi4e21tcsqa ...
            timestep: 0.0226
            t/T: 12.2
            breaking: 0
            Fi_iter: 4 Final_residual: 3.73e-09  Fi_time: 0.00174
            breaking: 0
            Fi_iter: 4 Final_residual: 3.08e-09  Fi_time: 0.00172
            breaking: 0
            Fi_iter: 5 Final_residual: 1.68e-11  Fi_time: 0.00167
            umax: 0.102
            vmax: 0
            wmax: 0.0941
            dt: 0.0226
            wavegentime: 2.49e-05
            reinitime: 0
            gctime: 0.000539         average gctime: 0.0005
            Xtime: 0.000126  average Xtime: 0.000641
            total time: 7.62451   average time: 0.00878
            timer per step: 0.00885
            ------------------------------------
            869
            simtime: 19.6
            timestep: 0.0226
            ...

            # Real-time streaming of a task's standard error output
            $ inductiva logs f0bnqgf4fcr4asgi4e21tcsqa --stderr
            
            # Wait for the task to start before streaming logs
            $ inductiva logs f0bnqgf4fcr4asgi4e21tcsqa --wait

            # Stream logs without colors
            $ inductiva logs f0bnqgf4fcr4asgi4e21tcsqa --no-color
    """)

    parser.set_defaults(func=stream_task_logs_tail)
