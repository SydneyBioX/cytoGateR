# Convert a tree to a tibble representation

Convert a tree to a tibble representation

## Usage

``` r
tree_to_df(node, parent = NA_character_, id = "root")
```

## Arguments

- node:

  A gating tree (node or leaf).

- parent:

  Parent node id (used internally for recursion).

- id:

  Node id for the current subtree (used internally for recursion).

## Value

A tibble with columns \`id\`, \`parent\`, \`type\`, \`label\`.

## Examples

``` r
tree <- list(
  type = "node", depth = 0, cells = 1:6, marker = "CD3e",
  cutoff = 0.5, sep_score = 0.8,
  left = list(type = "leaf", depth = 1, cells = 1:3),
  right = list(type = "leaf", depth = 1, cells = 4:6)
)
tree_to_df(tree)
#> # A tibble: 3 × 4
#>   id     parent type  label                  
#>   <chr>  <chr>  <chr> <chr>                  
#> 1 root   NA     node  "CD3e > 0.50\nsep=0.80"
#> 2 root_L root   leaf  "leaf\nn=3"            
#> 3 root_R root   leaf  "leaf\nn=3"            
```
