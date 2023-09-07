# Sample index file

This is a sample Markdown file that should be the index of the project
documentation, which should be built using [Sphinx](https://www.sphinx-doc.org).

Below are some examples on how to include common content types.

## Math formulas

Inline math formulas can be included using the `$...$` syntax from LaTeX.
For example, `$e^{i\pi} + 1 = 0$` is rendered as $e^{i\pi} + 1 = 0$.

Formulas spanning the full text width can be written using the `$$...$$` syntax.

For example,

```tex
$$
\sum_{k = 1}^n k  = \frac{n (n + 1)}{2}
$$
```

is rendered as

$$
\sum_{k = 1}^n k  = \frac{n (n + 1)}{2}
$$

## Code snippets

Common markdown syntax can be used to render code snippets:

```python
def square(x):
    """Computes the square of an input x, x^2."""
    return x*x
```

## Docstrings

Docstrings of modules, classes, functions, etc. can be automatically rendered in
a documentation file.

In **reStructuredText** (`.rst` files), docstrings can be included using the
following syntax:

```rst
.. automodule:: package.module
.. autoclass:: package.module.class
.. autofunction:: package.module.function
```

See the [documentation](https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html)
for a complete list of the available automatic docstring documentation
directives.

To embed docstrings in **Markdown** (`.md` or `.markdown` files), the [`eval-rst`
directive](https://myst-parser.readthedocs.io/en/latest/syntax/syntax.html#syntax-directives-parsing)
from [MyST](https://myst-parser.readthedocs.io/) should be used. For example, the directives in the
example above should be written as

````markdown
```{eval-rst}
.. automodule:: package.module
.. autoclass:: package.module.class
.. autofunction:: package.module.function
```
````

For more details, see this [usage guide](https://myst-parser.readthedocs.io/en/latest/sphinx/use.html#use-sphinx-ext-autodoc-in-markdown-files).
