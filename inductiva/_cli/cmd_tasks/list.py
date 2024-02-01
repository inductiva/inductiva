"""List the tasks information via CLI."""

from inductiva import tasks, utils


def list_tasks(args):
    """ List the last user's tasks. 

    Lists based on the flags (task_id, last_n).
    It can print the last N tasks (default value is 5)
    or a single task with a specific ID
    """
    if args.task_id is not None:
        task_list = [tasks.Task(args.task_id)]

    else:
        last_n = 5 if args.last_n is None else args.last_n
        task_list = tasks.get(last_n=last_n)

    table_dict = tasks.to_dict(task_list)

    def color_formater(status):
        if status == "success":
            return utils.format_utils.emphasis_formater(status, "green")
        elif status in ["failed", "killed", "executer-failed", "zombie"]:
            return utils.format_utils.emphasis_formater(status, "red")
        return status

    formatters = {
        "Submitted": [
            utils.format_utils.datetime_formatter,
            utils.format_utils.spacing_formater
        ],
        "Started": [
            utils.format_utils.datetime_formatter,
            utils.format_utils.spacing_formater
        ],
        "Status": [utils.format_utils.spacing_formater],
        "ID": [utils.format_utils.spacing_formater],
        "Simulator": [utils.format_utils.spacing_formater],
        "Computation Time": [utils.format_utils.spacing_formater],
    }

    header_formatters = [lambda x: x.upper()]
    if not utils.format_utils.getenv_bool("DISABLE_TERMINAL_EMPHASIS", False):
        formatters["Status"].insert(0, color_formater)
        header_formatters.append(
            lambda x: utils.format_utils.emphasis_formater(x, "bold"))

    print(
        utils.format_utils.get_tabular_str(table_dict,
                                           formatters=formatters,
                                           header_formatters=header_formatters))


def register(parser):
    """Register the list user's tasks command."""

    subparser = parser.add_parser("list", help="List tasks.")
    group = subparser.add_mutually_exclusive_group()
    group.add_argument("-n",
                       "--last-n",
                       type=int,
                       help="List last N tasks. Default: 5.")
    group.add_argument("-id",
                       "--task-id",
                       type=str,
                       help="List a task with a specific ID.")

    subparser.set_defaults(func=list_tasks)
