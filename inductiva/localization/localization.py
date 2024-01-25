"""Localization module that provides support for multiple languages."""
from typing import Optional
from threading import RLock
import pathlib
import logging
import locale
import json

_INDUCTIVA_DEFAULT_LANG = "en"

_logger = logging.getLogger(__name__)


class Translator:
    """Translation provider with string rendering support.

    Attributes:
        TRANSLATIONS_DIR (pathlib.Path): Path to the directory where the
            translation files are stored.
    """
    TRANSLATIONS_DIR = pathname = pathlib.Path(__file__).parent / "translations"

    class _Template(str):
        """Custom string class to allow for the use of the % operator"""

        def __mod__(self, other):
            """
            Allow for the use of the % operator to format the string
            with the str.format method.
            """
            if isinstance(other, dict):
                return self.format(**other)
            elif isinstance(other, tuple):
                return self.format(*other)
            return self.format(other)

        def __str__(self) -> str:  # pylint: disable=invalid-str-returned
            # hack to allow logging to work with this class.
            # logging uses msg = str(self.msg) in getMessage() before applying
            # the % operator. Here we return the object itself to avoid route
            # the string formatting to the __mod__ method
            return self

    def __init__(self, lang: Optional[str] = None) -> None:
        self._translations = {}
        self._lock = RLock()
        self._lang = None

        # make sure the default language is loaded
        # before setting it to the required value
        self.set_lang(_INDUCTIVA_DEFAULT_LANG)
        self.set_lang(lang)

    def get_lang(self):
        """Get the current language."""
        return self._lang

    def set_lang(self, lang: str):
        """Set the current language.

        Args:
            lang (str): Name of the language to set.
        """
        with self._lock:
            # if None is provided, use the system's locale language if
            # available. Otherwise use the packages's default.
            if lang is None:
                loc, _ = locale.getlocale()
                lang = loc.split("_")[0] if loc else _INDUCTIVA_DEFAULT_LANG

            if lang in self._translations:
                return

            self._load(lang)
            self._lang = lang

    def _load(self, lang: str) -> None:
        """Load the translation file for the given language.

        NOTE: This method is not thread-safe.

        Args:
            lang (str): Name of the language to load.
        """
        pathname = self.TRANSLATIONS_DIR / f"{lang}.json"

        if not pathname.exists():
            # if translation file not available, there is not point on trying
            # to use a translation for this message. Just use the default.
            _logger.debug(
                "Translation file not available for language '%s'. "
                "Will use '%s' as default language.", lang,
                _INDUCTIVA_DEFAULT_LANG)
            return

        with open(pathname, "r", encoding="utf-8") as f:
            translations = json.load(f)

        for key, translation in translations.items():
            if isinstance(translation, list):
                translation = "\n".join(translation)
            translations[key] = self._Template(translation)

        self._translations[lang] = translations

    def format(self, key: str, *args, **kwargs) -> str:
        """Return a formatted version of the translation for the given key
        in the current language, using substitutions from args and kwargs.

        If the key is not found in the currently set language, the default
        language will be used. If the key doesn't exist in the default language,
        the key itself will be used as the translation.

        Args:
            key (str): Key for the translation to format.
            *args: Positional arguments for the str.format method.
            **kwargs: Keyword arguments for the str.format method.
        """
        return self[key].format(*args, **kwargs)

    def __getitem__(self, key: str):
        """Get the translation for the given key in the currently language.

        If the key is not found in the current language, the default
        language will be used. If the key doesn't exist in the default language,
        the key itself will be returned.

        Args:
            key (str): Key for the translation to get.
        """

        translations = self._translations.get(self._lang)
        if key not in translations:
            translations = self._translations.get(_INDUCTIVA_DEFAULT_LANG)
        if key not in translations:
            return self._Template(key)
        return translations.get(key, key)


# pseudo-singleton object that is already initialized to the system's locale
# or the pakcage's default language if the system's locale is not available.
translator = Translator()
