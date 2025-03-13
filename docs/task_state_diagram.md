
```{mermaid}
---
title: Task state diagram
---
stateDiagram-v2
    direction TB
    PENDING_INPUT : PENDING_INPUT\n(Waiting for input upload) 
    ZOMBIE
    SUBMITTED
    KILLED

    [*] --> PENDING_INPUT : upon submission by user
    PENDING_INPUT --> SUBMITTED : after all input \n files are uploaded
    PENDING_INPUT --> KILLED : after being \n killed by the user
    PENDING_INPUT --> ZOMBIE : when machine \n group is terminated

    SUBMITTED --> STARTED : upon being picked \n by an executor
    SUBMITTED --> KILLED : after being killed  \n by the user
    SUBMITTED --> ZOMBIE : when machine group \n gets terminated
    ZOMBIE --> [*]


    STARTED --> PENDING_KILLED : user requested kill \n but is not yet delivered to executor
    PENDING_KILLED --> KILLED : when kill request \n is delivered to executor
    STARTED --> EX_TERM_USER : when machine group \n is terminated (by monitoring service or user)
    STARTED --> SPOT_PREEMPTED : when machine group \n is SPOT, and the executor gets preempted
    STARTED --> EX_FAILED : error in executor machine \n (eg no space in disk, etc)
    STARTED --> SUCCESS : Execution succeeded and \n output is available for download
    STARTED --> EX_TERM : when executor is \n terminated (not by user)
    STARTED --> FAILED : execution failed with \n non-success return code 
    
    EX_TERM --> SUBMITTED : when task is automatically \n resubmitted (hidden to user)
    SPOT_PREEMPTED --> SUBMITTED : when task is automatically \n resubmitted (hidden to user)
    

    EX_TERM_USER --> [*]
    EX_FAILED --> [*]
    FAILED --> [*]
    KILLED --> [*]
    SUCCESS --> [*]

```

