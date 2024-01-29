
def remove(args):
    """List storage contents."""
    path = args.path
    confirm = args.confirm
    if not confirm:
        prompt = input(f"Are you sure you want to remove {path}? (y/n)")
        confirm = prompt.lower() in ["y", "ye", "yes"]
    
    if confirm:
        print(f'removed {path}')

def register(parser):
    subparser = parser.add_parser("remove", aliases=["rm"],
        help="remove remote storage entries")
    subparser.add_argument("path", default="/")
    subparser.add_argument("-y", action="store_true", dest="confirm", default=False)
    subparser.set_defaults(func=remove)


