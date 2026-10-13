# Minimal spatial proteomics example

A small \`SpatialExperiment\` containing selected samples for
demonstrating cytoGateR visualization and annotation functions.

## Usage

``` r
data(cytoGateR_example)
```

## Format

A \`SpatialExperiment\` with 40 markers and the selected cells. The
object contains:

- assay:

  An \`"exprs"\` expression assay.

- colData:

  Sample identifiers, image names, reference labels, tree-gating labels
  and probabilities, high-confidence core labels, and weighted
  k-nearest-neighbour predictions, confidence values, and full
  \`KNN_P\_\*\` class-probability columns.

- spatialCoords:

  Two-dimensional cell coordinates.

## Source

A subset of the IMMUcan 2022 Cancer Example dataset from the
\`imcdatasets\` package.
