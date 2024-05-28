# Inductiva Documentation

## Deployment

The tutorials site can be deployed to staging or production from a source branch
via dedicated GitHub actions.

1. Commit and push your changes to a branch (can be a PR branch, main or
   something else)
2. Go over to the [Actions](https://github.com/inductiva/inductiva/actions) page
3. Click [Deploy Tutorials to Staging Server](https://github.com/inductiva/inductiva/actions/workflows/deploy-tutorials-to-staging.yml) or [Deploy Tutorials to Production Server](https://github.com/inductiva/inductiva/actions/workflows/deploy-tutorials-to-prod.yaml)
4. Click 'Run workflow ▼' button
5. Select the branch containing your changes
6. Click 'Run workflow'

## Documentation build tips

When working on the documentation it can be rather useful to run a local
build of the static site and have it reload when any changes are detected.

This can be performed in a very straightforward way using [sphinx-autobuild](https://github.com/sphinx-doc/sphinx-autobuild#readme)

### Set Python 3.8 as the local Python version for the build

We'll be using [PyEnv](https://github.com/pyenv/pyenv) for this first step but
feel free to use any other Python version manager.

```console
brew update
brew install pyenv
cd tutorials
pyenv install 3.8
pyenv local 3.8
```

### Install all sphinx dependencies including,`sphinx-autobuild`

```console
pip3 install -r requirements.txt
```

### Launch a local documentation build and monitor it for changes

```console
sphinx-autobuild . /tmp/inductiva-tutorials -W
```

* `sphinx-autobuild . /tmp/inductiva-tutorials` will build the tutorials and watch for
changes
* `-W` will handle [warnings as errors](https://www.sphinx-doc.org/en/master/man/sphinx-build.html#cmdoption-sphinx-build-W)
