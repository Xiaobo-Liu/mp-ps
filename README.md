# Mixed-Precision Paterson--Stockmeyer Method for Evaluating Polynomials of Matrices

Function `include/PS_mp.m` computes the Taylor approximants to the matrix exponential in arbitrary precision; function `include/cosm_hsd.m` is the counterpart that works in double precision and uses only double,
single, and half (simulated by `chop`) precisions (no `mp` computation). Function `include/PS_mp_gen.m` is for general polynomials of matrices and takes the scalar coefficients (of decaying modulus) as one of the inputs.

Details on the underlying algorithms can be found in the technical report:

X. Liu. Mixed-Precision Paterson--Stockmeyer Method for Evaluating Polynomials of Matrices, MIMS EPrint, Manchester Institute for Mathematical Sciences, The University of Manchester, UK, December 2023.

All codes used for generating the data in the above report are included in this repository.

## Dependencies

The code in this repository may require the Advanpix Multiprecision Computing Toolbox for MATLAB (www.advanpix.com).

## License

See `license.txt` for licensing information.
