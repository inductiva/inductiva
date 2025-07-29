"""Top command for task."""
import textwrap
from typing import TextIO
import argparse
import sys

from inductiva import tasks


def top(args: argparse.Namespace, fout: TextIO = sys.stdout):
    """Prints the result of the `top -b -H -n 1` command.
    
    This command will list the processes and threads (-H) in batch mode (-b).
    This command will run only once (-n 1) instead of running continuously.
    The result is an instant snapshot of the machine CPU and RAM metrics.
    """
    task_id = args.id
    task = tasks.Task(task_id)
    # pylint: disable=protected-access
    result, return_code = task._top()
    if result:
        print(result, file=fout)
    return return_code


def register(parser):
    """Register the `top` command for monitoring task resource usage."""
    subparser = parser.add_parser(
        "top",
        help="Displays the output of the `top` command from the machine running"
        " the task.",
        formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva tasks top` command streams the output of the "
        "`top -b -H -n 1` command executed on the machine where the "
        "task is running. This allows real-time monitoring of system "
        "processes and resource usage for the task.")

    subparser.add_argument("id",
                           type=str,
                           help="The ID of the task to monitor.")

    subparser.epilog = textwrap.dedent("""\
        examples:        
        $ inductiva tasks top qpusar8bch509k56g1hvv5yxk
        top - 12:00:15 up 18 min,  0 users,  load average: 1.14, 0.99, 0.58
        Threads: 226 total,   2 running, 224 sleeping,   0 stopped,   0 zombie
        %Cpu(s): 24.2 us,  1.5 sy,  0.0 ni, 72.7 id,  0.0 wa,  0.0 hi,  1.5 si,  0.0 st
        MiB Mem :  16008.2 total,  12976.4 free,   1057.1 used,   1974.7 buff/cache
        MiB Swap:      0.0 total,      0.0 free,      0.0 used.  14656.3 avail Mem 

            PID USER      PR  NI    VIRT    RES    SHR S  %CPU  %MEM     TIME+ COMMAND
            1469 task-ru+  20   0  894208 711108  36128 R  99.9   4.3   9:56.46 d_hydro+
            1557 task-ru+  20   0    9016   3812   3140 R   6.2   0.0   0:00.01 top
                1 root      20   0  165128  10828   7912 S   0.0   0.1   0:01.22 systemd
                2 root      20   0       0      0      0 S   0.0   0.0   0:00.00 kthreadd
                3 root       0 -20       0      0      0 I   0.0   0.0   0:00.00 rcu_gp
                4 root       0 -20       0      0      0 I   0.0   0.0   0:00.00 rcu_par+
        ...
    """)

    subparser.set_defaults(func=top)
