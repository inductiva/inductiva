# Access your Remote Storage

After your simulations finish, the results are securely stored in your **remote 
storage bucket**. The Inductiva CLI connects you to your remote bucket, and
provides you with various commands to effectively manage all the data of your 
simulation outputs.

## Check Storage Usage

You can check how much storage space your simulations are currently occupying:

```console
$ inductiva storage size
Total user's remote storage in use: 2.79 GB
```

## List Storage Contents

You can have a detailed view of what's in your storage, including sorting options 
for easier navigation:

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

## Clean Up Storage
Once you've backed up your data locally or need to free up space, you can delete 
data from your remote storage. You can delete several paths within your
storage, or even delete everything with the `remove` command.

Here's an example where you remove a path within your storage. Use this command 
with caution as it permanently deletes data from your remote 
storage: 

```console
$ inductiva storage remove hodbisrxjxhdbkknv60xmy6ti/
# The CLI always prompts for confirmation to prevent accidental data loss
You are about to remove the following paths from your remote storage space:
  - 0bet8jrpp2gz974n42nsd9n2p/
Are you sure you want to proceed (y/[N])? y
Removing '0bet8jrpp2gz974n42nsd9n2p/' from remote storage...
Successfully removed '0bet8jrpp2gz974n42nsd9n2p/' from remote storage.
```

You can alternatively clear all storage by adding the `--all` flag to the `remove` command.
Such commands will always be followed with a confirmation prompt to ensure user intention
and prevent irreversible loss.

These are the main functionalities of the CLI for storage at the moment. In case
you would like to have more functionalities via the CLI [contact us](mailto:support@inductiva.ai).
