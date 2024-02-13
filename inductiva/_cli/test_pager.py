"""
Example of using the PagedOutput class with a scheduled
update of the content.
"""
from ansi_pager import PagedOutput, StoppableScheduler
from datetime import datetime


def producer(fhandle: PagedOutput, schedulr: StoppableScheduler):
    fhandle.clear()
    count = schedulr.count
    lines = [f"--- line {count} {i} ---" for i in range(10 + count)]
    fhandle.writelines(lines)
    now = datetime.now().strftime("%H:%M:%S")
    fhandle.set_footer(f"iteration {count} at {now}")


every = 2.0  # seconds
header = f"> every {every} s: inductiva tasks list -n 10"

with PagedOutput(header, "footer") as fout:
    scheduler = StoppableScheduler(every, producer, args=())
    scheduler.args = (fout, scheduler)
    scheduler.start()
    fout.run()
    scheduler.stop()
    scheduler.join()
