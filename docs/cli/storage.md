# Remote Storage

Finally, when the simulation finishes the results are saved in the user's remote bucket. 

Hence, the CLI allows users to connect to their remote bucket where all the data of simulations live and manage it as they wish.

To start, they can explore the storage use at any time with:

```bash
$ inductiva storage size
Total user's remote storage in use: 0.46 GB
```

Thereafter, the contents of the storage can be listed as follows:
```bash
$ inductiva storage list --max-results 10 --order-by size --sort-order desc
Name                        Size      Creation Time
--------------------------  --------  ----------------
hodbisrxjxhdbkknv60xmy6ti/  87.55 MB  01 Feb, 14:55:43
m4d487kf46n4pp868qddzn67m/  54.48 MB  01 Feb, 16:49:36
ibeg639dv96yo2kx18wbtxsfi/  51.93 MB  01 Feb, 14:57:25
2oqv8tyfa1dubeq63z4g5e2zn/  36.67 MB  01 Feb, 16:48:38
4athv38xc7s17v79ucazzr0i0/  35.39 MB  01 Feb, 16:50:56
f0bnqgf4fcr4asgi4e21tcsqa/  12.32 MB  02 Feb, 00:58:28
opqo99i4axn6fde5pjvu5hzpt/  12.32 MB  02 Feb, 00:57:35
57mr4kas99jxb9titkeackano/  11.52 MB  01 Feb, 23:07:17
mak1ji62s7axf7mespkc36g7e/  11.52 MB  01 Feb, 23:07:14
ox8718m0pwfi02zczui3qky4w/  11.52 MB  01 Feb, 23:07:16
```

And when having downloaded the data to their local machines, or simply avoiding having too much clutter on the remote storage, users can quickly delete several paths within their storage or, if they wish, remove everything. Tread carefully with the following command, but in any case, you will be asked for confirmation:

```
$ inductiva storage remove hodbisrxjxhdbkknv60xmy6ti/
Are you sure you want to remove hodbisrxjxhdbkknv60xmy6ti/? (y/[N]) y
Removing hodbisrxjxhdbkknv60xmy6ti/ in the remote storage.
Successfully removed remote path 'hodbisrxjxhdbkknv60xmy6ti/'.
```

These are the main functionalities of the CLI at the moment. In case you would like to have more access via the CLI contact us at [contact].

#### What to read next
* []()
