"""List the storage contents via CLI."""

from inductiva import storage


def listdir(args):
    """List the user's remote storage contents."""
    storage.listdir(args.path, args.max_results, args.order_by, args.sort_order)


def register(parser):
    """Register the list user's remote storage command."""

    subparser = parser.add_parser("list",
                                  aliases=["ls"],
                                  help="List remote storage contents.")
    subparser.add_argument("path", default="/", type=str, nargs="?")
    subparser.add_argument("-m", "--max-results", default=10, type=int)
    subparser.add_argument("-o",
                           "--order-by",
                           default="creation_time",
                           type=str,
                           choices=["creation_time", "size"])
    subparser.add_argument("-s",
                           "--sort-order",
                           default="desc",
                           type=str,
                           choices=["desc", "asc"])

    subparser.set_defaults(func=listdir)
