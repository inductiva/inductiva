"""Compares a list of tasks via CLI."""
import argparse
import inductiva


def compare(args):
    """Kills a task by id."""
    ids = args.id
    sort_by = args.sort_by
    reverse = args.reverse
    no_info = 1 if args.no_info else 2
    ids = set(ids)

    list_of_tasks = [inductiva.tasks.Task(id) for id in ids]

    inductiva.tasks.compare(list_of_tasks, sort_by, reverse, verbose=no_info)

    return 0


def register(parser):
    """Register the compare tasks command."""

    subparser = parser.add_parser("compare",
                                  help="Compares a list of tasks.",
                                  formatter_class=argparse.RawTextHelpFormatter)

    subparser.description = (
        "The `inductiva tasks compare` command prints a comparison of tasks\n"
        "The comparison can be sorted by the table columns and the sort order\n"
        "can be reversed. In the end the information of the first task\n"
        "is printed.\n")

    subparser.add_argument("id",
                           type=str,
                           help="ID(s) of the task(s) to compare.",
                           nargs="+")
    subparser.add_argument("--sort-by",
                           type=str,
                           default="Computation",
                           help="Sort by a column.")
    subparser.add_argument("-r",
                           "--reverse",
                           action="store_true",
                           help="Determines the sort order.")
    subparser.add_argument(
        "--no-info",
        action="store_true",
        help="Do not print the information of the first task.")

    subparser.set_defaults(func=compare)
