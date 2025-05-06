# Step 3: Stop the Local Task-Runner
To stop the task-runner Docker container, press `Ctrl+C`.

Note that this does **not** terminate the associated machine group. To shut it down, use the Inductiva CLI:

```bash
$ inductiva resources terminate <machine-group-name>
```

Alternatively, you can terminate it using the Inductiva API:

```python
machine_group.terminate()
```
