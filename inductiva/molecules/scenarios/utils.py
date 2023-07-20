"""Post processing auxiliary functions."""
import time
import os
import io
import base64
import shutil
import threading
from PIL import Image
from moviepy.editor import ImageSequenceClip
from ipywidgets import Output, IntProgress


class MovieMaker:
    """ 
    Class to create a movie from a trajectory rendered on 
    NGLview.
    """

    def __init__(
        self,
        view,
        output_path='movie.mp4',
        fps=8,
        start=0,
        stop=-1,
        step=1,
        timeout=0.1,
    ):
        """
        Initializes a MovieMaker object to construct a movie from the frames 
        rendered on NGLview. 
        Args:
            view: NGLview object. 
            output_path: Path to save the movie.
            fps: Frame rate (in frames per second).
            start: Starting frame number.
            stop: Ending frame number. 
            step: Step size for the frames.
            timeout: The waiting time between rendering two consecutive frames.
        """

        self.view = view
        self.fps = fps
        self.timeout = timeout
        self.output_path = output_path

        if stop < 0:
            stop = self.view.max_frame + 1

        self._range = range(start, stop, step)
        self._event = threading.Event()
        self.frames_dir = 'frames'
        os.makedirs('frames', exist_ok=True)

    def sleep(self):
        time.sleep(self.timeout)

    def make(self):
        self.progress = IntProgress(description='Rendering...',
                                    max=len(self._range) - 1)
        self._event = threading.Event()

        def _make(event):
            image_files = []
            iw = None
            for i in self._range:
                self.progress.value = i
                if not event.is_set():
                    self.view.frame = i
                    self.sleep()
                    iw = self.view.render_image()
                    self.sleep()
                    im_bytes = base64.b64decode(self.view._image_data)
                    im_bytes = io.BytesIO(im_bytes)
                    try:
                        image = Image.open(im_bytes)
                        filename = os.path.join(self.frames_dir,
                                                f'figure_{i}.png').format(i)
                        image_files.append(filename)
                        image.save(filename, "PNG")
                    except Exception:
                        print("Could not save frame", i)
                        continue

                    if iw:
                        iw.close()

            if not self._event.is_set():
                self.progress.description = "Writing ..."
                clip = ImageSequenceClip(image_files, fps=self.fps)
                with Output():
                    if self.output_path.endswith('.gif'):
                        clip.write_gif(
                            self.output_path,
                            fps=self.fps,
                            verbose=False,
                        )
                    else:
                        clip.write_videofile(self.output_path, fps=self.fps)
                self.progress.description = 'Done'
                time.sleep(1)
                self.progress.close()

        self.thread = threading.Thread(target=_make, args=(self._event,))
        self.thread.daemon = True
        self.thread.start()
        shutil.rmtree(self.frames_dir)
        return self.progress
