"""Localization module that provides support for multiple languages."""
from typing import Optional
from threading import RLock
import pathlib
import logging
import locale
import json

_INDUCTIVA_DEFAULT_LANG = "en"

_logger = logging.getLogger(__name__)


class Template:
    """Template for lazy string formatting.

    Formatting is performed using the str.format method only when a string
    representation of the object is required (typically using str(obj)).
    """

    def __init__(self, fmt: str, /, *args, **kwargs):
        """Intialize a Template object using the given formatting string
        and the positional and keyword substituion arguments.

        Args:
            fmt (str): The template string.
            *args: Positional arguments for the str.format method.
            **kwargs: Keyword arguments for the str.format method.
        """
        self.fmt = fmt
        self.args = args
        self.kwargs = kwargs

    def __str__(self) -> str:
        """Return a formatted version of self.fmt, using substitutions
        from self.args and self.kwargs."""
        return self.fmt.format(*self.args, **self.kwargs)


class Translator:
    """Translation provider with string rendering support.

    Attributes:
        TRANSLATIONS_DIR (pathlib.Path): Path to the directory where the
            translation files are stored.
    """

    TRANSLATIONS_DIR = pathname = pathlib.Path(__file__).parent / "translations"

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

            if self._load(lang):
                self._lang = lang

    def _load(self, lang: str) -> bool:
        """Load the translation file for the given language.

        Returns True if the language file was successfully
        loaded, False otherwise.
        NOTE: This method is not thread-safe.

        Args:
            lang (str): Name of the language to load.
        """
        if not lang:
            return False

        lang = lang.lower()
        pathname = self.TRANSLATIONS_DIR / f"{lang}.json"

        if not pathname.exists():
            # if translation file not available, there is not point in trying
            # to use a translation for this message. Just use the default.
            _logger.debug(
                "Translation file not available for language '%s'. "
                "Will use '%s' as default language.", lang,
                _INDUCTIVA_DEFAULT_LANG)
            return False

        with open(pathname, "r", encoding="utf-8") as f:
            translations = json.load(f)

        for key, translation in translations.items():
            if isinstance(translation, list):
                translations[key] = "\n".join(translation)

        self._translations[lang] = translations
        return True

    def __call__(self, key, *args, **kwargs) -> Template:
        """Return a Template object with the translation for the given key as
        formatting string, and substitutions from args and kwargs.

        The translation is retrieved from the active language. If the key is
        not found in the active language, the default language will be used.
        If the key doesn't exist in the default language, the key itself will
        be used as the translation.

        Args:
            key (str): Key for the translation which will be used as formatting
                string in the Template constructor.
            *args: Positional substitutions for the Template constructor.
            **kwargs: Keyword substitutions for the Template constructor.
        """
        translation = self[key]
        return Template(translation, *args, **kwargs)

    def __getitem__(self, key: str) -> str:
        """Get the translation for the given key in the active language.

        If the key is not found in the active language, the default language
        will be used. If the key doesn't exist in the default language, the key
        itself will be returned.

        Args:
            key (str): Key for the translation to get.
        """
        translations = self._translations.get(self._lang)
        if translations is None or key not in translations:
            translations = self._translations.get(_INDUCTIVA_DEFAULT_LANG)
        return translations.get(key, key)


# pseudo-singleton object that is already initialized to the system's locale
# or the package's default language if the system's locale is not available.
translator = Translator()
