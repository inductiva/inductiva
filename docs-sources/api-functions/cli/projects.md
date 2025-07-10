# projects

The `inductiva projects` command provides utilities to manage and
retrieve information about your projects. You can list all your projects
and download associated task files.

## Usage

```sh
inductiva projects [-h] {download,list} ...
```

### Description
Projects in Inductiva serve as folders for organizing tasks.
Using the CLI, you can list your projects and retrieve the
outputs of all tasks in a specific projects.

## Options

- **`-h, --help`** → Show help message and exit.

## Available Subcommands

### `list (ls)`
List all existing projects associated with your account.

```sh
inductiva projects list
```

or using the alias:

```sh
inductiva projects ls
```

### `download`
Download the task files of a specific project.

```sh
inductiva projects download <project_id>
```

## Example Usage

### List all projects:
```sh
inductiva projects list
```

### Download the task files for a specific project:
```sh
inductiva projects download my_project_id
```

## Need Help?
Run the following command for more details:

```sh
inductiva projects --help
```
