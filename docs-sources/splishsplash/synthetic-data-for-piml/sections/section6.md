# Visualizing Multiple Simulations in a Single Animation

Now that we generated multiple animations, it's time to bring the results together into a single visualization. This will help us clearly compare how different parameters affect particle movement over time.

We’ll do this by combining individual simulation GIFs into a single animated grid. The result is one comprehensive animation that showcases multiple simulations side-by-side.

## Creating a Combined GIF

The following Python script takes several GIF files and arranges them in a grid layout to produce a single, unified GIF:

```python
from PIL import Image, ImageSequence
import math
import argparse

def combine_gifs_grid(gif_paths, output_path, columns, rows=None, max_seconds=None, bg_color=(255, 255, 255, 0), frame_duration=100):
    num_gifs = len(gif_paths)

    # Determine number of rows automatically if not specified
    if rows is None:
        rows = math.ceil(num_gifs / columns)
    elif rows * columns < num_gifs:
        print(f"Warning: grid size {rows}x{columns} is too small for {num_gifs} GIFs. Some will be skipped.")

    # Load GIF frames
    gifs_frames = []
    gifs_sizes = []
    max_frames = 0

    for path in gif_paths:
        im = Image.open(path)
        frames = [frame.copy().convert("RGBA") for frame in ImageSequence.Iterator(im)]
        gifs_frames.append(frames)
        gifs_sizes.append(im.size)
        max_frames = max(max_frames, len(frames))

    # Adjust frame count based on time limit
    if max_seconds is not None:
        max_frames = min(max_frames, int((max_seconds * 1000) / frame_duration))

    max_width = max(w for w, h in gifs_sizes)
    max_height = max(h for w, h in gifs_sizes)
    out_width = max_width * columns
    out_height = max_height * rows

    combined_frames = []

    for frame_idx in range(max_frames):
        new_frame = Image.new("RGBA", (out_width, out_height), bg_color)

        for i, frames in enumerate(gifs_frames):
            if i >= rows * columns:
                break

            frame = frames[frame_idx % len(frames)]
            x = (i % columns) * max_width
            y = (i // columns) * max_height
            new_frame.paste(frame, (x, y), frame)

        combined_frames.append(new_frame)

    combined_frames[0].save(
        output_path,
        save_all=True,
        append_images=combined_frames[1:],
        duration=frame_duration,
        loop=0,
        disposal=2,
        transparency=0,
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine multiple GIFs into a grid layout.")
    parser.add_argument("output", help="Filename of the output GIF")
    parser.add_argument("columns", type=int, help="Number of columns in the grid")
    parser.add_argument("rows_or_gifs", nargs='+', help="Optionally specify the number of rows, followed by the input GIF paths")
    parser.add_argument("--seconds", type=float, help="Limit the output GIF duration (in seconds)")

    args = parser.parse_args()

    try:
        rows = int(args.rows_or_gifs[0])
        gif_paths = args.rows_or_gifs[1:]
    except ValueError:
        rows = None
        gif_paths = args.rows_or_gifs

    combine_gifs_grid(gif_paths, args.output, args.columns, rows, args.seconds)
    print(f"Combined GIF saved as {args.output}")
```

> 💡 **Note**: Make sure to install the required dependencies by running:
>
> ```bash
> pip install pillow
> ```

## Example Usage

You can run the script with the following command:

```bash
python combine_gifs.py combined.gif 3 3 gifs/output_0.gif gifs/output_1.gif gifs/output_2.gif gifs/output_3.gif gifs/output_4.gif gifs/output_5.gif gifs/output_6.gif gifs/output_7.gif gifs/output_8.gif --seconds 7
```

### What this command does:

* `combined.gif`: Name of the output file.
* `3`: Number of columns in the GIF grid.
* `3`: Number of rows in the GIF grid.
* `gifs/output_*.gif`: Input GIFs from previous simulations.
* `--seconds 7`: (Optional) Limit the total duration of the animation to 7 seconds.

This will generate a 3×3 animation grid from 9 individual simulation outputs, giving you a unified view of how different configurations perform.

<p align="center"><img src="../../_static/combined.gif" alt="Visualization of 9 simulations" width="600"></p>

---

Ready to take your simulation visualizations to the next level? Let's go!
