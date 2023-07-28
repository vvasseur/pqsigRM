[Enhanced pqsigRM](https://csrc.nist.gov/Projects/pqc-dig-sig/round-1-additional-signatures) is a hash and sign scheme which belongs to the family of GPV-like signatures but using error correcting codes instead of lattices.
In this framework the independence of the signatures distribution from the secret key is critical for security

We discovered significant biases in the signature distribution of enhanced pqsigRM, which expose information about the secret key.


Thomas Debris-Alazard, Pierre Loisel and Valentin Vasseur


# Bias

We provide an initial script that:
* generates a key,
* produces 1k signatures,
* identifies matched pairs,
* counts the number of recovered matched pairs (using the secret key).

This script shows that signatures are biased.

You can execute this by running `correlation.sh`

## Execution time

Under a minute.

## Dependencies

The same as those of enhanced pqsigRM reference implementation, with the addition of Python.


# Recovering the (U|U+V) structure of the secret key

We have a second script that:
* generates a key,
* produces 100k signatures,
* identifies matched pairs,
* reveals the (U|U+V) structure of the code.

You can execute this by running `attack.sh`

## Execution time

Approximately 3 hours.

## Dependencies

The same as those of enhanced pqsigRM reference implementation, and SageMath is also required.
