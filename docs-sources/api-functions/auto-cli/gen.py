from pathlib import Path

LEVELS_TO_STR = {
    1: '-',
    2: '~',
    3: '^',
}


def extract_cmd(cmd_dir: Path) -> str:
    return cmd_dir.name.removeprefix("cmd_").replace('_', '-')

def get_cmd_dirs(root: Path) -> str:
    return sorted([
        d for d in root.iterdir() 
        if d.is_dir() and d.name.startswith("cmd_")
    ])


def gen_cmd_block(header: str, command: str, level: int, file):
    underline = LEVELS_TO_STR[level] * len(header)
    file.write(f"""
{header}
{underline}

.. program-output:: inductiva {command} --help
   :prompt:
""")


def gen_cmd_blocks(cmd_dir: Path, cmd: str, level: int, file):
    sub_cmds = sorted([
        f.stem.replace('_', '-') for f in cmd_dir.iterdir() 
        if f.is_file() and
            f.stem != "__init__" and
            f.name.endswith(".py")
    ])

    for sub_cmd in sub_cmds:
        gen_cmd_block(f"``{sub_cmd}``", f"{cmd} {sub_cmd}", level, file)

    cmd_dirs = get_cmd_dirs(cmd_dir)
    for cmd_dir in cmd_dirs:
        cmd_name = extract_cmd(cmd_dir)
        cmd += " " + cmd_name
        gen_cmd_block(f"``{cmd_name}``", cmd, level, file)
        gen_cmd_blocks(cmd_dir, cmd, level + 1, file)


def gen_file(cmd_dir: Path):
    cmd = extract_cmd(cmd_dir)
    filename = f"{cmd}.rst"
    with open(filename, "w", encoding="utf-8") as file:
        header = f"inductiva {cmd}"
        underline = '=' * len(header)
        file.write(f"{header}\n{underline}\n")
        gen_cmd_block("Overview", cmd, 1, file)
        gen_cmd_blocks(cmd_dir, cmd, 1, file)


def gen_files(cli_dir: str):
    cli_path = Path(cli_dir)
    cmd_dirs = get_cmd_dirs(cli_path)
    for cmd_dir in cmd_dirs:
        gen_file(cmd_dir)


CLI_DIRECTORY = "../../../inductiva/_cli"
gen_files(CLI_DIRECTORY)
