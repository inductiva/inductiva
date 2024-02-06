from typing import Optional, Type, Iterator, AnyStr, Iterable
from types import TracebackType
import curses.ascii
import curses
import threading
import typing
import time
import io


class PagedOutput(io.TextIOBase):
    SCROLL_DOWN = 1
    SCROLL_HOME = 0
    SCROLL_UP = -1

    HEADER = 1
    FOOTER = 1

    def __init__(self, name=None):
        self._lock = threading.RLock()
        self._init()
        self._set_header(name)
        self.flush()

    def _init(self):
        """Intialize the display."""
        self._window = curses.initscr()
        self._current_line = 0
        self._line_count = 0
        self._height, self._width = self._window.getmaxyx()

        curses.noecho()
        curses.cbreak()
        # hide the cursor
        curses.curs_set(0)

        # receive keypad events and make 1-millisecond
        # timeouts on `getch`
        self._window.keypad(True)
        self._window.timeout(1)

        self._txtpad_view_height = self._height - self.HEADER - self.FOOTER
        self._header = curses.newpad(self.HEADER, self._width)
        self._footer = curses.newpad(self.FOOTER, self._width)
        self._txtpad = curses.newpad(1000, self._width)

        self._header.bkgd(" ", curses.A_BOLD)
        self._footer.bkgd(" ", curses.A_BOLD)

        self._update_footer()
        self.flush()

    def _set_header(self, content: str):
        """Set the header of the display.

        This method is typically called with the name of the file as content.
        Note: The display is only updated after a call to the `flush` method.
        """
        self._header.erase()
        if content:
            self._header.addstr(content)

    def _update_footer(self):
        """Update the line tracking information on the footer pad."""
        with self._lock:
            curr_line = self._current_line
            if self._line_count:
                curr_line += 1

            left = "  q or ESC to quit, ← ↑ → ↓ to navigate"
            right = f"{curr_line}/{self._line_count}  "
            fill = " " * (self._width - len(left) - len(right) - 1)

            self._footer.erase()
            self._footer.addstr(left + fill + right)

    def _scroll(self, direction, page=False):
        """Scroll the content of the text pad in the given direction.

        This method scrolls a single line in the given direction when
        `page` is False. Otherwise, an entire page is scrolled.
        """
        with self._lock:
            line = self._current_line
            if page:
                line += direction * self._txtpad_view_height
            else:
                line += direction
            self.scroll_to(line)

    def scroll_to(self, line):
        """Scroll the display to the given line.

        The method scrolls the display so that the given line is at the top
        of the text pad. The line is clamped to the range [0, line_count-1].
        """
        with self._lock:
            line = min(self._line_count - 1, max(0, line))
            self._current_line = line
            self._update_footer()

    def clear(self) -> None:
        """Clear the content of the file-like object.

        Both the text and footer pads are cleared and the display refreshed.
        """
        with self._lock:
            self._current_line = 0
            self._line_count = 0

            self._update_footer()
            self._txtpad.erase()

            self.flush()

    def run(self):
        """Run the navigation loop"""
        try:
            while True:
                self.flush()
                ch = self._window.getch()

                if ch in (curses.ascii.ESC, 113):
                    # q or ESC
                    break
                elif ch == curses.KEY_HOME:
                    self._scroll(self.SCROLL_HOME, page=True)
                elif ch == curses.KEY_DOWN:
                    self._scroll(self.SCROLL_DOWN)
                elif ch == curses.KEY_UP:
                    self._scroll(self.SCROLL_UP)
                elif ch in (curses.KEY_LEFT, 98):
                    # left or b
                    self._scroll(self.SCROLL_UP, page=True)
                elif ch in (curses.KEY_RIGHT, curses.ascii.SP):
                    # right or space
                    self._scroll(self.SCROLL_DOWN, page=True)

        except KeyboardInterrupt:
            return

    def write(self, text: str) -> int:
        """Write string to the file-like object.

        Returns the total number of lines in the display.
        Note: The display is only updated after a call to the `flush` method.
        """
        with self._lock:
            self._txtpad.addstr(text)

            self._line_count += text.count("\n")
            current_line = self._line_count - self._txtpad_view_height
            self._current_line = max(current_line, 0)
            self._update_footer()
            return self._line_count

    def writelines(self, lines: Iterable[AnyStr]) -> None:
        """Write a list of lines to the file-like object.

        Line separators are not added, so it is usual for each of the
        lines provided to have a line separator at the end.
        Note: The display is only updated after a call to the `flush` method.
        """
        self.write("".join(lines))

    def flush(self) -> None:
        """Flush the content of the file-like object to the display."""
        with self._lock:
            self._header.refresh(0, 0, 0, 0, self.HEADER, self._width)

            self._txtpad.refresh(self._current_line, 0, self.HEADER, 0,
                                 self._height - self.FOOTER - 1, self._width)

            self._footer.refresh(0, 0, self._height - self.FOOTER, 0,
                                 self._height, self._width)

    def close(self) -> None:
        """Close the file-like object.

        The method closes the file-like object and resets the terminal
        to its previous state."""
        self._window.keypad(False)
        curses.nocbreak()
        curses.echo()
        curses.endwin()

    def __enter__(self) -> "PagedOutput":
        return self

    def __exit__(self, t: Optional[Type[BaseException]],
                 value: Optional[BaseException],
                 traceback: Optional[TracebackType]) -> Optional[bool]:
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
        # self._txtpad.move(offset, 0)
        # self._current_line = offset
        raise io.UnsupportedOperation("File-like object is not seekable")

    def tell(self) -> int:
        raise io.UnsupportedOperation(
            "File-like object does not support the tell method")

    def truncate(self, size: Optional[int] = ...) -> int:
        # self._txtpad.clrtobot()
        # self._line_count = self._current_line
        raise io.UnsupportedOperation(
            "File-like object does not support the truncate method")

    def __next__(self) -> AnyStr:
        raise TypeError("File-like object is not an iterator")

    def __iter__(self) -> Iterator[AnyStr]:
        raise TypeError("File-like object is not iterable")


class Scheduler(threading.Thread):

    def __init__(self, every: float, action, args=(), kwargs=None):
        super().__init__()
        self.every = every
        self.action = action
        self.args = args
        self.kwargs = kwargs or {}

        self._event = threading.Event()
        self._count = 0

    @property
    def count(self):
        """Return the number of times the action has been executed."""
        return self._count

    def run(self):
        """Run the scheduler until the stop method is called."""
        elapsed = 0.0
        while not self._event.is_set():
            self._count += 1

            now = time.time()
            self.action(*self.args, **self.kwargs)
            elapsed = time.time() - now

            # take into account the time the action takes to complete
            if elapsed < self.every:
                self._event.wait(self.every - elapsed)

    def stop(self):
        """Stop the scheduler."""
        self._event.set()


def main():
    every = 5.0  # seconds

    def producer(fout: PagedOutput, schedulr: Scheduler):
        fout.clear()
        for i in range(50):
            fout.write(f"--- line {schedulr.count} {i+1} ---\n")
        fout.scroll_to(0)

    with PagedOutput(f"> every {every} s: inductiva tasks list -n 10") as f:
        scheduler = Scheduler(every, producer, args=())
        scheduler.args = (f, scheduler)
        scheduler.start()
        f.run()
        scheduler.stop()
        scheduler.join()


if __name__ == "__main__":
    main()
