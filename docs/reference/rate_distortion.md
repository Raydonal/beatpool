# Rate–distortion curve of a quantised spectral codec

Quantises the coefficients of a linear projection with a uniform scalar
quantiser whose step is a multiple of the coefficient standard
deviation, entropy-codes the integer levels, and reports the mean
relative distortion and the rate in bits per beat. The basis and the
mean are charged once per data holder and amortised over `n_amortize`
beats, which is the cost model used in the paper.

## Usage

``` r
rate_distortion(
  X,
  T,
  k = ncol(T),
  steps = 2^seq(0, -8, by = -0.5),
  center = TRUE,
  sample_bits = 11L,
  basis_bits = 16L,
  mean_bits = 16L,
  n_amortize = nrow(X),
  transmit_basis = TRUE
)
```

## Arguments

- X:

  A beat matrix, one beat per row.

- T:

  A basis matrix.

- k:

  Number of components.

- steps:

  Multipliers of the coefficient standard deviation; a larger multiplier
  gives a coarser quantiser, a lower rate and a larger distortion.

- center:

  Center before projecting.

- sample_bits:

  Bits per raw sample in the reference codec, used only for the
  compression ratio.

- basis_bits:

  Bits per entry of the transmitted basis.

- mean_bits:

  Bits per entry of the transmitted mean.

- n_amortize:

  Number of beats over which the basis and the mean are amortised;
  defaults to the number of beats in `X`.

- transmit_basis:

  Whether the basis is transmitted. Fixed transforms such as the DCT are
  known to the decoder and cost nothing.

## Value

A data frame with the step multiplier `step`, the rate in bits per beat,
the mean relative distortion in percent and the compression ratio
against the raw reference.
