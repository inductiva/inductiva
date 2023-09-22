"""Tasks CLI subcommands."""
import inductiva


def list_tasks(args):
    inductiva.tasks.list(last_n=args.last_n)


def register_tasks_cli(parser):
    subparsers = parser.add_subparsers()
    list_subparser = subparsers.add_parser("list", help="List tasks.")
    list_subparser.add_argument(
        "-n",
        "--last-n",
        type=int,
        default=5,
        help="List last N tasks. Default: %(default)s",
    )
    # Register function to call when this subcommand is used
    list_subparser.set_defaults(func=list_tasks)
