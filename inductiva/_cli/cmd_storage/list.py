
def list(args):
    print(f'storage2 list with args {args}')

def register(parser):
    subparser = parser.add_parser("list", aliases=["ls"],
        help="List currently active resources")
    subparser.add_argument("path", default="/", nargs='?')
    subparser.add_argument("-o", "--order-by", choices=["size", "creation_time"])
    subparser.add_argument("-s", "--sort-order", choices=["asc", "desc"])

    subparser.set_defaults(func=list)


