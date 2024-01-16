"""Tasks CLI subcommands."""
from inductiva import tasks, utils
from inductiva import _cli


def list_tasks(args):
    """List tasks."""

    override_col_space = {
        "Submitted": 18,
        "Started": 18,
        "Status": 10,
        "Computation Time": 20,
        "Resource Type": 20
    }
    print(
        utils.format_utils.get_dataframe_str(
            tasks.list(last_n=args.last_n),
            override_col_space=override_col_space))


def register_tasks_cli(parser):
    _cli.utils.show_help_msg(parser)

    subparsers = parser.add_subparsers()

    list_subparser = subparsers.add_parser("list", help="List tasks")
    list_subparser.add_argument(
        "-n",
        "--last-n",
        type=int,
        default=5,
        help="List last N tasks. Default: %(default)s",
    )
    # Register function to call when this subcommand is used
    list_subparser.set_defaults(func=list_tasks)
