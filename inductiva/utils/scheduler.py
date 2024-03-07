"""Utility to run a function at regular time intervals."""
from typing import Union, Callable
import threading
import time


class StoppableScheduler(threading.Thread):
    """A stoppable scheduler for running a function at regular intervals.
    
    The scheduler runs a given function at regular intervals until the stop
    method is called. If the call to the action takes longer than the
    interval, the scheduler will wait for the remaining time before calling
    the action at the next interval.
    """

    def __init__(self,
                 every: Union[int, float],
                 action: Callable[..., None],
                 args=(),
                 kwargs=None):
        """Initialize the scheduler.

        Args:
            every (float): The interval at which the action should be executed.
            action (callable): The function to be executed.
            args (tuple, optional): The positional arguments to be passed to the
                action. Defaults to ().
            kwargs (dict, optional): The keyword arguments to be passed to the
                action. Defaults to None.
        """
        super().__init__()
        self.every = every
        self.action = action
        self.args = args
        self.kwargs = kwargs or {}

        self._event = threading.Event()
        self._count = 0

    @property
    def count(self) -> int:
        """Return the number of times the action has been executed."""
        return self._count

    def run(self) -> None:
        """Run the scheduler until the stop method is called."""
        every = self.every

        while not self._event.is_set():
            now = time.time()
            self._count += 1
            self.action(*self.args, **self.kwargs)
            elapsed = time.time() - now

            # take into account the time the action takes to complete
            #if elapsed < self.every:
            #    self._event.wait(self.every - elapsed)
            self._event.wait(every - elapsed % every)

    def stop(self) -> None:
        """Stop the scheduler."""
        self._event.set()
