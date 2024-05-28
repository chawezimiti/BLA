## Resubmission

This is a resubmission. In this version I have:

* Removed the extra spaces in the Description field of the DESCRIPTION file.

* Changed how the information message is writen to console for the functions
  R/exp_boundary.R; R/ble_profiling.R; R/StartValues.R; R/na_drop.R. I have now used the
  message() instead of cat() function.

* Used the on.exit() function to avoid changing user's options, par or working directory
  when they run the functions R/exp_boundary.R; R/summastat.R.


## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

- This is a first submission and hence can not eliminate the new release note.

