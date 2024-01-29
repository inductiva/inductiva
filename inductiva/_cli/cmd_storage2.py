from inductiva._cli import utils

def list(args):
    print(f'storage2 list with args {args}')

def remove(args):
    """List storage contents."""
    path = args.path
    confirm = args.confirm
    if not confirm:
        prompt = input(f"Are you sure you want to remove {path}? (y/n)")
        confirm = prompt.lower() in ["y", "ye", "yes"]
    
    if confirm:
        print(f'removed {path}')

def register(root_parser):

    parser = root_parser.add_parser('storage2',
                help='storage2 management commands')
    utils.show_help_msg(parser)

    parser = parser.add_subparsers()

    # LIST -------------
    subparser = parser.add_parser("list", aliases=["ls"],
        help="List currently active resources")
    subparser.add_argument("path", default="/", nargs='?')
    subparser.add_argument("-o", "--order-by", choices=["size", "creation_time"])
    subparser.add_argument("-s", "--sort-order", choices=["asc", "desc"])
    subparser.set_defaults(func=list)



    # REMOVE -------------
    subparser = parser.add_parser("remove", aliases=["rm"],
        help="remove remote storage entries")
    subparser.add_argument("path", default="/")
    subparser.add_argument("-y", action="store_true", dest="confirm", default=False)
    subparser.set_defaults(func=remove)
