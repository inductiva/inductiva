# How to Receive an Email When a Task Finishes  

## The Challenge  

Running simulations in the cloud means you don't have to wait around watching
progress bars, but how do you know when a task is done?

Without notifications, it's easy to lose track of task completions, leading to
delays in your workflow or missed opportunities to launch follow-up jobs.

## The Solution  

Inductiva lets you register an **event trigger** that sends you an email as soon
as a task finishes and its outputs are uploaded.  

You’ll receive the email **regardless of whether the task succeeded, failed, or timed out**, and the message will include the final task status, so you know exactly what happened.  

Here’s how to set it up:

```python
from inductiva import events

events.register(
    trigger=events.triggers.TaskOutputUploaded(task_id="<task_id>"),
    action=events.actions.EmailNotification(email_address="<your_email>")
)
```

Just replace `<task_id>` and `<your_email>` with the task id and your email address.  

Now, every time that task finishes, you’ll be notified immediately, no more refreshing dashboards or guessing when it’s done.  

Want to track multiple tasks? Simply register a notification for each one.