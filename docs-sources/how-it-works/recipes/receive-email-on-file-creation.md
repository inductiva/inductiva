# How to Receive an Email When a File Appears in a Simulation  

## The Challenge  

Some simulations don’t fail in an obvious way, there’s no red error message or
crash. Instead, they quietly generate a file like `error.log` or bury issues
inside a long output file. If you’re not actively checking for these signs,
you could miss problems entirely.  

Other times, you may simply want to track progress by waiting for a specific
file to be generated, like an intermediate result, status flag, or checkpoint.

## The Solution  

With Inductiva’s **observer triggers**, you can automatically detect when a
specific file is created **or** when a file contains a particular string (using
regex), and get notified by email right away.  

This allows you to monitor internal simulation events, catch subtle errors, or
just track key milestones—without manually digging through logs.

Here’s how to do it:

```python
from inductiva import events

# Notify when a specific file is created
events.register(
    trigger=events.triggers.ObserverFileExists(
        task_id=task.id,
        file_path="error.log"),
    action=events.actions.EmailNotification(
        email_address="your@email.com")
)

# Notify when a file contains a specific string (e.g. error, warning, or success pattern)
events.register(
    trigger=events.triggers.ObserverFileRegex(
        task_id=task.id,
        file_path="output.log",
        regex="ERROR"),
    action=events.actions.EmailNotification(
        email_address="your@email.com")
)
```

Replace `"your@email.com"` with your actual address, and customize the file path and
regex as needed.  

## Use Cases  

- 🔥 **Error detection**: Some simulators create an error file instead of failing outright.
- 📄 **Log monitoring**: Catch keywords like `ERROR` and `WARNING` inside output logs.
- ⏱️ **Progress tracking**: Get notified when key simulation stages finish by watching for intermediate files.
- ✅ **Completion signal**: Some simulations might create a "done" flag file—watch for that to know when to post-process.

Observer triggers let you stay one step ahead without constantly checking in.