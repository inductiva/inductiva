# FAQ

#### How can I know the overall size and cost of my remote storage?
You can check how much storage space you are currently using by issuing:

```sh
inductiva storage size
```
The result will inform you of both the size and the cost per/month of your storage:
```sh
Total storage size used:
	Volume: 41.98 GB
	Cost: 0.78 US$/month
```


#### How do I delete **ALL** my storage permanently?
To permanently delete all your storage you can add the `--all` flag to the `remove` command.
```sh
inductiva remove --all
```
This commands is be followed by a confirmation prompt to ensure the user intention
and prevent irreversible loss of data. All data will be permanently delete.


#### How to to continuously monitor task progress using the CLI?
You can do it by continuously calling `inductiva tasks list` for a specific task. For that, you can employ a 
[watch method](https://www.geeksforgeeks.org/watch-command-in-linux-with-examples/) on Linux or Mac, 
which will invoke the command `inductiva tasks list`  at set intervals and refresh task information.

In th example below, we refreshe task information every 10 seconds:

```sh
# Monitor task status updates every 10 seconds
$ watch -n 10 inductiva tasks list --id jxwt0rm8s8xspdfcegtgkkana
Every 10.0s: inductiva tasks list --id jxwt0rm8s8xspdfcegtgkkana                                                                                 


       ID                              SIMULATOR          STATUS         SUBMITTED              STARTED                COMPUTATION TIME         RESOURCE TYPE
       jxwt0rm8s8xspdfcegtgkkana       splishsplash       success        08 Feb, 13:25:49       08 Feb, 13:26:04       0:00:35                  c2-standard-4
```


#### Can I save my task logs locally?
You can save the log stream content and the stream consumer's status for later 
inspection by redirecting them to a file:

- **Redirect stdout to a file**: You can redirect stdout (file descriptor 1) to 
a file with the redirection operator (`>`) to save the entire log stream content. 
For example:

    ```bash
    $ inductiva logs TASK_ID 1>out.txt
    ```
- **Redirect stderr to a file**: You can redirect the status of the stream consumer, 
outputted to stderr (file descriptor 2), to a file like this:

    ```bash
    $ inductiva logs TASK_ID 2>err.txt
    ```
- **Redirect stdout and stderr to separate files**: You can redirect both stdout 
and stderr simultaneously to separate files using:

    ```bash
    $ inductiva logs TASK_ID 1>out.txt 2>err.txt
    ```
- **Disable ANSI globally**: If you need to disable ANSI escape codes, which are 
used to support the status bar at the bottom of the log stream, either export the 
`ANSI_ENABLED` environment variable or set it locally:

    ```bash
    $ export ANSI_ENABLED=0 # Applies to the entire shell session
    $ ANSI_ENABLED=0 inductiva logs TASK_ID. # Applies only to this command
    ```

However, you don't need to do so. You can visit the Task's page on the Web Console,
and you will find the logs of the task there ready to be downloaded.
