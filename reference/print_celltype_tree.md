# Print a cell-type gating tree

Print a cell-type gating tree

## Usage

``` r
print_celltype_tree(node, indent = "")
```

## Arguments

- node:

  A gating tree (node or leaf) produced by
  \[build_fullcoverage_tree()\].

- indent:

  String used internally for indentation during recursive printing.

## Value

Invisibly returns \`NULL\`.

## Examples

``` r
tree <- list(
  type = "node", depth = 0, cells = 1:6, marker = "CD3e",
  cutoff = 0.5, sep_score = 0.8,
  left = list(type = "leaf", depth = 1, cells = 1:3),
  right = list(type = "leaf", depth = 1, cells = 4:6)
)
print_celltype_tree(tree)
#> [depth 0] CD3e > 0.500  (sep = 0.80, n = 6)
#>  |- low
#>  |  [leaf] depth =1 n_cells =3
#>  `- high
#>     [leaf] depth =1 n_cells =3
```
