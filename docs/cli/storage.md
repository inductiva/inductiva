# Remote Storage

Finally, when the simulation finishes the results are saved in the user's remote bucket. 

Hence, the CLI allows users to connect to their remote bucket where all the data of simulations live and manage it as they wish.

**Explore size:**

To start, they can explore the storage use at any time with:

```bash
$ inductiva storage size
Total user's remote storage in use: 2.79 GB
```

**List contents:**

Thereafter, the contents of the storage can be listed as follows:
```bash
$ inductiva storage list --max-results 10 --order-by size --sort-order desc

       NAME                             SIZE           CREATION TIME
       0bet8jrpp2gz974n42nsd9n2p/       56.11 MB       06 Feb, 11:32:29
       05ujj5m0ytdkckxwk1tq1b5io/       27.93 MB       08 Feb, 09:19:44
       6a2h1wnxywpea8jfoxiikdjf7/       26.49 MB       07 Feb, 13:47:03
       f8joznwc9xf9a4nypcaei6v2s/       12.79 MB       07 Feb, 09:16:55
       dpq2cv6b5f9p1c77nc8anjo10/       12.00 MB       08 Feb, 09:39:31
       r4kerxf4b53krgn0s3fyece3b/       11.92 MB       07 Feb, 11:47:48
       j9qzrpiohgt7x97od3tw4wccd/       11.74 MB       07 Feb, 11:47:46
       iqi71gonoacfj7fknox3rvnq2/       11.52 MB       07 Feb, 11:47:45
       dxmnxdrfrv84pfbzbvm9v0dat/       11.43 MB       07 Feb, 11:47:43
       bgtwgnnyq5qa5hecegzdx6okr/       11.36 MB       07 Feb, 11:47:40
```

**Remove contents:**

When having downloaded the data to their local machines, or simply avoiding
having too much clutter on the remote storage, users can quickly delete several
paths within their storage or, if they wish, remove everything. Tread carefully
with the following command, but in any case, you will be asked for confirmation:

```
$ inductiva storage remove hodbisrxjxhdbkknv60xmy6ti/
You are about to remove the following paths from your remote storage space:
  - 0bet8jrpp2gz974n42nsd9n2p/
Are you sure you want to proceed (y/[N])? y
Removing 0bet8jrpp2gz974n42nsd9n2p/ in the user's remote storage.
Successfully removed remote path '0bet8jrpp2gz974n42nsd9n2p/'.
```

Notice that the `remove` command can take multiple paths as arguments, and even
remove everything with the `--all` flag.

These are the main functionalities of the CLI for storage at the moment. In case
you would like to have more functionalities via the CLI [contact us](mailto:support@inductiva.ai).

#### What to read next
* []()
