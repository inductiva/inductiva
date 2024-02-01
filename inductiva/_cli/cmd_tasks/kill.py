"""Kills a tasks by id via CLI."""

import inductiva


def kill_task(args):
    """Kills a task by id."""
    if not args.yes:
        prompt = input("Do you confirm you want to "
                       f"kill {len(args.id)} tasks (y/[N])?")
        confirm = prompt.lower() in ["y", "ye", "yes"]
        if confirm:
            print("Killing tasks.")
        else:
            print("Aborting the killing of tasks.")
            return

    for task_id in args.id:
        print("Killing task ", task_id)
        inductiva.tasks.Task(task_id).kill(wait_timeout=args.wait_timeout)


def register(parser):
    """Register the kill task command."""

    subparser = parser.add_parser("kill", help="Kill a running task.")
    subparser.add_argument("id",
                           type=str,
                           help="ID(s) of the task to kill.",
                           nargs='+')
    subparser.add_argument("-w",
                           "--wait-timeout",
                           type=float,
                           default=None,
                           help="Timeout to wait for the kill to complete.")
    subparser.add_argument("-y",
                           "--yes",
                           action='store_true',
                           help="Skip kill confirmation.")

    subparser.set_defaults(func=kill_task)
