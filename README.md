### Fortran codes for the field-field correlator of the McGuire model at finite temperatures

Firstly, choose the size and the range of the quadrature: 
```bash
gfortran legendre_rule.f90 && ./a.out
```

Then set the systems parameters within the main program and run it:
```bash
gfortran rho_x_exact.f90 && ./a.out
```

And similarly procedures apply to other files. 

For additional details, see [arXiv:2202.07657](https://arxiv.org/abs/2202.07657) and [arXiv:2308.06482](https://arxiv.org/abs/2308.06482).
