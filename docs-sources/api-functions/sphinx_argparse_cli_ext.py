"""Custom Sphinx extension to generate CLI documentation automatically."""

# pylint: disable=protected-access

import argparse
import re
import textwrap
from typing import Optional
from docutils import nodes
from sphinx_argparse_cli._logic import SphinxArgparseCli

from inductiva._cli.main import get_main_parser
from banner_small_directive import BannerSmallDirective


def _get_subparsers(
    parser: argparse.ArgumentParser,) -> Optional[argparse._SubParsersAction]:
    return next(
        (action for action in parser._actions
         if isinstance(action, argparse._SubParsersAction)),
        None,
    )


def _get_subparser(
    parser: argparse.ArgumentParser,
    name: str,
) -> argparse.ArgumentParser:
    subparsers = _get_subparsers(parser)
    return subparsers.choices[name]


def get_parser(command: str) -> argparse.ArgumentParser:
    subcommands = command.split(" ")
    parser = get_main_parser()
    for subcommand in subcommands:
        parser = _get_subparser(parser, subcommand)
    return parser


def get_parser_wrapper(command):

    def wrapper():
        return get_parser(command)

    return wrapper


class SphinxArgParseCliExt(SphinxArgparseCli):
    """Custom Sphinx extension to generate CLI documentation automatically."""

    RE_JOIN_LINES = re.compile(r"(?<!\n)\n(?!\n)(?=[A-Za-z])")
    RE_INLINE_CODE = re.compile(r"`([^`]+)`")

    def __init__(
        self,
        name: str,
        arguments: list[str],
        options: dict[str, Optional[str]],
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

        import sphinx_argparse_cli_ext  # pylint: disable=import-outside-toplevel
        setattr(sphinx_argparse_cli_ext, func_name, wrapper)

        options["module"] = __name__
        options["func"] = func_name
        options["usage_first"] = None
        options["epilog"] = ""

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

        origin = f"cli-{command}"
        self.banner = BannerSmallDirective.html(origin)

    def insert_transitions(self, root: nodes.Node):
        for section in root.findall(SphinxArgParseCliExt.is_options_section):
            parent = section.parent
            index = parent.index(section)
            parent.insert(index + 1, nodes.transition())

    @staticmethod
    def is_options_section(node: nodes.Node) -> bool:
        return (isinstance(node, nodes.section) and
                (title := node.next_node(nodes.title)) and
                "options" in title.astext().lower())

    @staticmethod
    def is_arguments_section(node: nodes.Node) -> bool:
        return (isinstance(node, nodes.section) and
                (title := node.next_node(nodes.title)) and
                "positional arguments" in title.astext().lower())

    @staticmethod
    def create_section(parser: argparse.ArgumentParser) -> nodes.Node:
        text = parser.epilog
        if text.startswith("examples:"):
            text = text.removeprefix("examples:")
        text = textwrap.dedent(text)

        code_block = nodes.literal_block(text, text, language="bash")

        reference = parser.prog.replace(" ", "-") + "-examples"
        section = nodes.section(ids=[reference])
        section += nodes.title(text="Examples")
        section += code_block

        return section

    def insert_examples_sections(self, root: nodes.Node):
        for section in root.findall(SphinxArgParseCliExt.is_options_section):
            title = section.next_node(nodes.title).astext()
            command = title.removeprefix("inductiva ").removesuffix(" options")

            parser = get_parser(command)
            if not parser.epilog:
                continue

            new_section = SphinxArgParseCliExt.create_section(parser)
            parent = section.parent
            index = parent.index(section)
            parent.insert(index + 1, new_section)

    def shorten_titles(self, root: nodes.Node):

        def _shorten_titles(condition, new_title):
            for section in root.findall(condition):
                title_node = section.next_node(nodes.title)
                new_title_node = nodes.title(text=new_title)
                title_node.parent.replace(title_node, new_title_node)

        _shorten_titles(
            SphinxArgParseCliExt.is_arguments_section,
            "Positional Arguments",
        )

        _shorten_titles(
            SphinxArgParseCliExt.is_options_section,
            "Options",
        )

    @staticmethod
    def format_paragraph(text_node: nodes.Text):
        text = text_node.astext()
        new_text = SphinxArgParseCliExt.RE_JOIN_LINES.sub(" ", text)
        text_node.parent.replace(text_node, nodes.Text(new_text))

    @staticmethod
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

        self.insert_transitions(root)
        self.insert_examples_sections(root)
        self.shorten_titles(root)

        formatters = [
            SphinxArgParseCliExt.format_paragraph,
            SphinxArgParseCliExt.format_inline_code,
        ]
        for formatter in formatters:
            self.format_text(root, formatter)

        return nodes_list + [self.banner]


def setup(app):
    app.add_directive("sphinx_argparse_cli_ext", SphinxArgParseCliExt)
    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
