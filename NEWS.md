# Changes in cit version 2.3.0 (2021-05-19)

* Fixed memory leaks with gsl (freeing memory).

* Fixed a bug in cit.bp() when using C. It was creating design matrices of
  the wrong size. That was causing memory leaks.

* Reamp for returning to CRAN.

* Adding a few tests for cit.bp().

