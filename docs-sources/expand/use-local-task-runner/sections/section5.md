# FAQ

## 1. Why do I get a permission error when running the task runner on Ubuntu 24.04?
Ubuntu 24.04 has a stricter AppArmor policy that restricts unprivileged user namespaces by default, which may prevent the task-runner from working.

To fix this, run the following command:

```
sudo sysctl -w kernel.apparmor_restrict_unprivileged_userns=0
```

## 2. What's the recommended Ubuntu version for running the task-runner?
We recommend Ubuntu 22.04 LTS for better compatibility and stability. If using Ubuntu 24.04, apply the fix above.
