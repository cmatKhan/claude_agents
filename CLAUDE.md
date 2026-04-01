# Global Preferences

## Output Style

Do not use emoji, arrows, bullets composed of special Unicode characters, or
other decorative non-ASCII characters in any output, including responses,
comments, and code. Plain punctuation (periods, commas, colons, dashes, etc.)
is fine.

## Documentation and Comments

Comments and docstrings should be informative but brief, and written in a
formal tone. Avoid redundant restatements of what the code plainly shows.

### Python

Use Sphinx-style docstrings for all functions and classes.

```python
def example(x, y):
    """Compute the sum of two values.

    :param x: First operand.
    :type x: float
    :param y: Second operand.
    :type y: float
    :returns: Sum of x and y.
    :rtype: float
    """
```

### R

Use roxygen2 documentation for all functions.

```r
#' Compute the sum of two values.
#'
#' @param x First operand (numeric).
#' @param y Second operand (numeric).
#' @return Numeric sum of x and y.
#' @export
example <- function(x, y) {
  x + y
}
```
