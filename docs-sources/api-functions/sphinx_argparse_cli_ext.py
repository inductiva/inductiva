import argparse
import re
import textwrap
from typing import Optional
from docutils import nodes
from sphinx_argparse_cli._logic import SphinxArgparseCli

from inductiva._cli.main import get_main_parser


def get_subparsers(
    parser: argparse.ArgumentParser,) -> Optional[argparse._SubParsersAction]:
    return next(
        (action for action in parser._actions
         if isinstance(action, argparse._SubParsersAction)),
        None,
    )


def get_subparser(
    parser: argparse.ArgumentParser,
    name: str,
) -> argparse.ArgumentParser:
    subparsers = get_subparsers(parser)
    return subparsers.choices[name]


def get_parser(command: str) -> argparse.ArgumentParser:
    subcommands = command.split(" ")
    parser = get_main_parser()
    for subcommand in subcommands:
        parser = get_subparser(parser, subcommand)
    return parser


def get_parser_wrapper(command):

    def wrapper():
        return get_parser(command)

    return wrapper


class SphinxArgParseCliExt(SphinxArgparseCli):

    RE_JOIN_LINES = re.compile(r'(?<!\n)\n(?!\n)(?=[A-Za-z])')
    RE_INLINE_CODE = re.compile(r'`([^`]+)`')

    def __init__(
        self,
        name: str,
        arguments: list[str],
        options: dict[str, str | None],
        content,
        lineno: int,
        content_offset: int,
        block_text: str,
        state,
        state_machine,
    ):
        command = content.pop()
        func_name = command.replace(" ", "_") + "_parser"
        wrapper = get_parser_wrapper(command)

        import sphinx_argparse_cli_ext
        setattr(sphinx_argparse_cli_ext, func_name, wrapper)

        options["module"] = __name__
        options["func"] = func_name

        super().__init__(
            name,
            arguments,
            options,
            content,
            lineno,
            content_offset,
            block_text,
            state,
            state_machine,
        )

    def is_options_section(node: nodes.Node) -> bool:
        return (isinstance(node, nodes.section) and
                (title := node.next_node(nodes.title)) and
                "options" in title.astext().lower())

    def create_section(text: str) -> nodes.Node:
        if text.startswith("examples:"):
            text = text.removeprefix("examples:")
        text = textwrap.dedent(text)

        code_block = nodes.literal_block(text, text, language="bash")

        section = nodes.section(ids=['examples-section'])
        section += nodes.title(text="Examples")
        section += code_block

        return section

    def insert_examples_sections(self, root: nodes.Node):
        subparsers = get_subparsers(self.parser)
        parsers = {str(self.parser): self.parser, **subparsers.choices}
        parsers = {n: p for n, p in parsers.items() if len(n) > 2}
        sections = root.findall(SphinxArgParseCliExt.is_options_section)

        for subcommand, node in zip(parsers, sections):
            parser = parsers[subcommand]
            if not parser.epilog:
                continue

            section = SphinxArgParseCliExt.create_section(parser.epilog)
            parent = node.parent
            index = parent.index(node)
            parent.insert(index + 1, section)

    def format_paragraph(text_node: nodes.Text):
        text = text_node.astext()
        new_text = SphinxArgParseCliExt.RE_JOIN_LINES.sub(' ', text)
        text_node.parent.replace(text_node, nodes.Text(new_text))

    def format_inline_code(text_node: nodes.Text):
        text = text_node.astext()
        parts = SphinxArgParseCliExt.RE_INLINE_CODE.split(text)
        if len(parts) == 1:
            return
        new_children = [
            nodes.Text(part) if i % 2 == 0 else nodes.literal(text=part)
            for i, part in enumerate(parts)
        ]
        text_node.parent.replace(text_node, new_children)

    def format_text(self, root: nodes.Node, formatter):
        for paragraph_node in root.findall(nodes.paragraph):
            for text_node in paragraph_node.findall(nodes.Text):
                formatter(text_node)

    def run(self):
        nodes_list = super().run()
        assert len(nodes_list) == 1
        root = nodes_list[0]

        self.insert_examples_sections(root)

        formatters = [
            SphinxArgParseCliExt.format_paragraph,
            SphinxArgParseCliExt.format_inline_code,
        ]
        for formatter in formatters:
            self.format_text(root, formatter)

        return nodes_list


def setup(app):
    app.add_directive("sphinx_argparse_cli_ext", SphinxArgParseCliExt)
    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
