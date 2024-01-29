import importlib
import pkgutil
import logging

_logger = logging.getLogger(__name__)

def load_commands(parser, module_dir, package, prefix=''):
    modules = pkgutil.iter_modules([module_dir])
    for _, name, _ in modules:

        if not name.startswith(prefix):
            continue

        modname = f'.{name}'
        
        try:
            module = importlib.import_module(modname, package)
        except Exception as e:
            _logger.error(f"Failed to import {modname}: {e}")
            continue
        
        if hasattr(module, 'register'):
            n = package.count('.')
            indent = n * '  '
            print(f"{indent}Registering {modname}")
            _logger.debug(f"%sRegistering %s",
                          package.count('.') * ' ', modname)
            module.register(parser)

