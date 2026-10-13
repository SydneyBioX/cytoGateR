# Define a small feedforward classifier

Returns raw logits (not softmax) – torch's cross-entropy loss expects
logits.

## Usage

``` r
.build_dl_net(input_dim, hidden_dims, n_classes, dropout)
```
