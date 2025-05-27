"""A simple ANSI-based pager for displaying formatted text in a terminal.

The module provides a file-like object that can be used to display text in a
terminal. The text is displayed in a paged manner, allowing the user to scroll
up and down, and to quit the display.
The pager requires the terminal to support ANSI escape codes.
"""
from typing import Optional, Type, Iterator, AnyStr, Iterable
from collections import namedtuple
from types import TracebackType
from functools import partial
import threading
import termios
import typing
import shutil
import tty
import sys
import re
import io
import os

from inductiva.utils import scheduler

strip_escaped = partial(re.compile(r"\x1b\[[;\d]*[A-Za-z]").sub, "")
strip_escaped.__doc__ = "Return the string with all ANSI escape codes removed."


def get_output_size():
    """Return the size of the terminal window as a named tuple."""
    return shutil.get_terminal_size()


Pad = namedtuple("Pad", "width height column line")
Pad.__doc__ = """A named tuple representing a pad in the terminal.
the tuple has the following fields:
- width: The width of the pad.
- height: The height of the pad.
- column: The column where the pad starts.
- line: The line where the pad starts.
"""


class Display:
    """A display object that holds the pads for the header, body and footer."""
    __slots__ = ["window", "header", "body", "footer"]

    def __init__(self):
        """Initialize the display object."""
        self.set_size(*get_output_size())

    def set_size(self, width: int, height: int) -> None:
        """Set the size of the display.
        Sets the size of the window and updates the sizes of the
        header, body and footer pads.

        Args:
            width (int): The width of the display.
            height (int): The height of the display.
        """
        self.window = Pad(width, height, 1, 1)
        self.header = Pad(width, 1, 1, 1)
        self.body = Pad(width, height - 2, 1, 2)
        self.footer = Pad(width, 1, 1, height)

    def resized(self) -> bool:
        """Return whether the terminal has been resized.

        If the terminal has been resized, the size of the display is updated and
        all the pads are recalculated. The method returns True if the terminal
        has been resized, False otherwise.
        """
        size = get_output_size()
        resized = size != (self.window.width, self.window.height)
        if resized:
            self.set_size(*size)
        return resized


class PagedOutput(io.TextIOBase):
    """A file-like object for displaying text in a terminal.

    The object provides a simple interface for displaying text in a terminal
    in a paged manner. The text is displayed in a buffer and the user can
    navigate through the buffer using the arrow keys. The object is also
    capable of displaying a header and a footer, both of which can be configured
    by the user to display custom content.
    """
    DEFAULT_FOOTER = "Type q to quit, ← ↑ → ↓ to navigate."

    SCROLL_DOWN = 1
    SCROLL_HOME = 0
    SCROLL_UP = -1

    def __init__(self, header: str = None, footer: str = DEFAULT_FOOTER):
        self._lock = threading.RLock()
        self._display = Display()
        self._scrolled_to = 0
        self._resizer = None
        self._header = None
        self._footer = None
        self._buffer = []

        self._enter_alt_screen()
        self.set_header(header)
        self.set_footer(footer)

        self._resizer = scheduler.StoppableScheduler(0.5, self._check_resize)
        self._resizer.daemon = True
        self._resizer.start()

    def _check_resize(self) -> None:
        """Check if the terminal has been resized and update the display."""
        with self._lock:
            if self._display.resized():
                self._render_header()
                self._render_body(self._scrolled_to)
                self._render_footer()

    def _enter_alt_screen(self) -> None:
        """Enter alternate screen and hide cursor."""
        sys.stdout.write("\033[?1049h\033[?25l")

    def _exit_alt_screen(self) -> None:
        """Exit alternate screen and show cursor."""
        sys.stdout.write("\033[?1049l\033[?25h")

    def _moveto(self, pad: Pad) -> None:
        """Move cursor to the start of the given pad."""
        sys.stdout.write(f"\033[{pad.line};{pad.column}H")

    def _clear_line(self) -> None:
        """Clear line where cursor is located."""
        sys.stdout.write("\033[2K")

    def _clear_screen(self) -> None:
        """Clear screen from cursor position till the end."""
        sys.stdout.write("\033[0J")

    def _render_header(self) -> None:
        """Render the header onto the display.

        Renders the content of the _header attribute to the header line
        onto the display.
        """
        self._moveto(self._display.header)
        self._clear_line()
        if self._header:
            sys.stdout.write(self._header)
        sys.stdout.flush()

    def _render_body(self, start: int) -> None:
        """Render the content of the text pad onto the display.

        Renders the content of the buffer, starting from the given line.
        """
        with self._lock:
            if start > self._scrolled_to + self._display.body.height:
                return

            self._moveto(self._display.body)
            self._clear_screen()
            start = self._scrolled_to
            end = start + self._display.body.height
            for line in self._buffer[start:end]:
                sys.stdout.write(line + "\n")

            sys.stdout.flush()

    def _render_footer(self) -> None:
        """Render the footer onto the display.

        Renders the content of the _footer + line tracking information
        to the footer line on the display
        """
        line_count = len(self._buffer)
        curr_line = self._scrolled_to
        if line_count:
            curr_line += 1

        left = self._footer if self._footer else ""
        stripped_left = strip_escaped(left)
        right = f"{curr_line}/{line_count}"
        fill = self._display.footer.width - len(stripped_left) - len(right)

        self._moveto(self._display.footer)
        self._clear_line()
        sys.stdout.write(left + fill * " " + right)
        sys.stdout.flush()

    def _scroll(self, direction: int, page: bool = False) -> None:
        """Scroll the content of the text pad in the given direction.

        This method scrolls a single line in the given direction when
        `page` is False. Otherwise, an entire page is scrolled.
        """
        with self._lock:
            line = self._scrolled_to
            if page:
                line += direction * self._display.body.height
            else:
                line += direction
            self.scroll_to(line)

    def set_header(self, header: str) -> None:
        """Set the header of the display.

        Args:
            header (str): The header to be displayed.
        """
        with self._lock:
            self._header = header
            self._render_header()

    def set_footer(self, footer: str) -> None:
        """Set the footer of the display.

        Args:
            footer (str): The footer to be displayed."""
        with self._lock:
            self._footer = footer if footer is not None else self.DEFAULT_FOOTER
            self._render_footer()

    def scroll_to(self, line: int) -> None:
        """Scroll the display to the given line.

        The method scrolls the display so that the given line is at the top
        of the text pad. The line is clamped to the range [0, line_count-1].
        """
        with self._lock:
            scrolled_to = min(len(self._buffer) - 1, max(0, line))
            self._scrolled_to = max(0, scrolled_to)
            self._render_body(self._scrolled_to)
            self._render_footer()

    def clear(self) -> None:
        """Clear the content of the file-like object.

        Both the text and footer pads are cleared and the display refreshed.
        """
        with self._lock:
            self._scrolled_to = 0
            self._buffer.clear()
            self._moveto(self._display.body)
            self._clear_screen()
            self._render_footer()

    def run(self) -> None:
        """Run the navigation loop"""

        stdin = sys.stdin.fileno()
        tattr = termios.tcgetattr(stdin)
        tty.setcbreak(stdin, termios.TCSANOW)

        scroll = {
            "\033[A": partial(self._scroll, self.SCROLL_UP, False),
            "\033[B": partial(self._scroll, self.SCROLL_DOWN, False),
            "\033[C": partial(self._scroll, self.SCROLL_DOWN, True),
            "\033[D": partial(self._scroll, self.SCROLL_UP, True),
            " ": partial(self._scroll, self.SCROLL_DOWN, True),
            "b": partial(self._scroll, self.SCROLL_UP, True),
        }

        try:
            buf = code = ""
            while True:
                buf = sys.stdin.read(1)
                if buf == "q":
                    break
                if action := scroll.get(buf):
                    action()
                elif buf == "\033":
                    code = buf
                elif buf in "[ABCD":
                    code += buf
                else:
                    code = ""
                    continue

                if action := scroll.get(code):
                    action()
                    code = ""

        except KeyboardInterrupt:
            termios.tcsetattr(stdin, termios.TCSANOW, tattr)

    def write(self, text: str) -> int:
        """Write string to the file-like object.

        Returns the total number of lines in the display.
        Note: The display is only updated after a call to the `flush` method.
        """
        with self._lock:
            start = len(self._buffer)
            lines = text.splitlines()
            self._buffer.extend(lines)
            self._render_header()
            self._render_body(start)
            self._render_footer()
            return len(self._buffer) - start

    def writelines(self, lines: Iterable[AnyStr]) -> None:
        """Write a list of lines to the file-like object.

        Line separators are not added, so it is usual for each of the
        lines provided to have a line separator at the end.
        Note: The display is only updated after a call to the `flush` method.
        """
        with self._lock:
            start = len(self._buffer)
            self._buffer.extend(x for line in lines for x in line.splitlines())
            self._render_body(start)
            self._render_footer()

    def flush(self) -> None:
        """Flush the content of the file-like object to the display."""
        sys.stdout.flush()

    def close(self) -> None:
        """Close the file-like object.

        The method closes the file-like object and resets the terminal
        to its previous state."""
        self._resizer.stop()
        self._resizer.join()
        self._exit_alt_screen()
        os.system("stty sane")  # Required to return terminal to normal.

    def __enter__(self) -> "PagedOutput":
        return self

    def __exit__(self, unused_t: Optional[Type[BaseException]],
                 unused_value: Optional[BaseException],
                 unused_traceback: Optional[TracebackType]) -> Optional[bool]:
        self.close()

    def isatty(self) -> bool:
        """Return whether this is an 'interactive' stream.

        Always returns False.
        """
        return False

    def writable(self) -> bool:
        """Return whether the object can be written to.

        Always returns True.
        """
        return True

    def readable(self) -> bool:
        """Return whether the object can be read from.

        Always returns False.
        """
        return False

    def seekable(self) -> bool:
        """Return whether the object supports random access.

        Always returns False.
        """
        return False

    def fileno(self) -> int:
        raise io.UnsupportedOperation("Underlying file object has no fileno")

    def read(self, n: int = ...) -> AnyStr:
        raise io.UnsupportedOperation("Underlying file object is not readable")

    def readline(self, limit: int = ...) -> AnyStr:
        raise io.UnsupportedOperation("Underlying file object is not readable")

    def readlines(self, hint: int = ...) -> typing.List[AnyStr]:
        raise io.UnsupportedOperation("Underlying file object is not readable")

    def seek(self, offset: int, whence: int = ...) -> int:
        raise io.UnsupportedOperation("File-like object is not seekable")

    def tell(self) -> int:
        raise io.UnsupportedOperation(
            "File-like object does not support the tell method")

    def truncate(self, size: Optional[int] = ...) -> int:
        raise io.UnsupportedOperation(
            "File-like object does not support the truncate method")

    def __next__(self) -> AnyStr:
        raise TypeError("File-like object is not an iterator")

    def __iter__(self) -> Iterator[AnyStr]:
        raise TypeError("File-like object is not iterable")
