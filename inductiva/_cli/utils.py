"""CLI utils."""


def show_help_msg(parser):
    """Show help message for command if no subcommand is provided."""
    parser.set_defaults(func=lambda _: parser.print_help())
