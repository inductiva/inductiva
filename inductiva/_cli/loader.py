"""Loader module to register all modules within each package automatically."""
import importlib
import pkgutil
import logging

_logger = logging.getLogger(__name__)


def load_commands(parser,
                  module_dir,
                  package,
                  prefix="",
                  ignores_prefix=None,
                  hides_prefix=None):
    """Load all commands within the module_dir and package."""
    modules = pkgutil.iter_modules([module_dir])
    metavar = []
    for _, name, _ in modules:
        if not name.startswith(prefix) or \
            (ignores_prefix is not None and name.startswith(ignores_prefix)):
            continue

        modname = f".{name}"
        try:
            module = importlib.import_module(modname, package)
        except Exception as e:  # pylint: disable=broad-except
            _logger.error("Failed to import %s: %s", modname, e)
            continue

        if hasattr(module, "register"):
            indent = package.count(".") * "  "
            _logger.debug("%sRegistering %s", indent, modname)
            module.register(parser)
            metavar.append(name)

    if hasattr(parser, "choices") and hides_prefix:
        parser.metavar = (
            "{" +
            ",".join(cmd for cmd, name in zip(parser._name_parser_map, metavar)  # pylint: disable=protected-access
                     if not name.startswith(hides_prefix)) + "}")
