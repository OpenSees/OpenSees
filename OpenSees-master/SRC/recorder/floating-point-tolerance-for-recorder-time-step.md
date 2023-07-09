# Selection of floating point tolerance for recorder time step

There are two precision errors that need to be balanced. 
The two options for setting a tolerance are evaluated below with code snippets and theoretical limits, 
and the simple scalar multiplier of the recorder time step has less error if chosen correctly. 
The number could be a user input but there is little scope that the number could have, maybe within 1e-4 to 1e-8, 
but the user must have very good control of the precision error that is developing between recordings, 
whereas 1e-5 is a quite reasonable default value. I have provided a full explanation of the different errors below. 
I have also made a suggestion for an addition to the documentation to explain the range when the recorder is reliable. 
There is also a third important precision error that is related to analyses in general, not to recorders.

Since the line is called for every recorder for every analysis step, it has been kept it as a single 
line with minimal calculation steps.

## Precision errors

tolerance check: `timeStamp - nextTimeStampToRecord >= -tol`

### First precision error: 
`timeStamp` incurs a maximum error every analysis step of the size `err = timeStamp * DBL_EPSILON`, 
where `DBL_EPSILON` is the smallest difference between 1 and a representable number that is not 1, 
and typically has a magnitude close to 1e-15 although can be larger or smaller depending on 
the operating system. The error in `timeStamp` can be cumulative 
(worst case scenario - and often is if time step is constant) therefore `n` analysis steps after 
the last recording the maximum error between `timeStamp` and `nextTimeStampToRecord` is
 `n * timeStamp * DBL_EPSILON`. If the analysis time step, `adt`, is very small then `n` will 
 be large and the error will be large. If the error is larger than the defined floating point 
 tolerance, `tol`, then an additional step will be taken before recording, which will mean the 
 error in the time recording will be the size of `adt - tol`, which can be as big as `deltaT`, 
 thus potentially 100% error (this often happens in the current implementation since there is no 
 floating point tolerance).

### Second precision error: 
This error is when the analysis time step is smaller than `tol`, and therefore records earlier than 
it should. The potential error on this is the size of `tol`.

### Third precision error: 
`timeStamp` incurs a maximum error every analysis step of the size `err = timeStamp * DBL_EPSILON` this 
can be cumulative over the whole analysis and result in a time drift. Therefore a limit 
of `adt > 1 /g * timeStamp * DBL_EPSILON * n_analysis` should be used, where `g` is the tolerable 
error in total time divided by the analysis time step, which ideally would be less than 1% by maybe 
be as high as 100% when `n_analysis` and `timeStamp` are large. Assuming that the user wishes to limit 
the drift to 10% over one million steps, then `adt > 1/.1 * timeStamp * DB_EPSILON * 1e6`. 
Therefore `adt > timeStamp * DB_EPSILON * 1e7`. Essentially very small time steps cannot be used 
when the current time is very large.


## Tolerance options:

Since the error is governed by both the size of `timeStamp` and `deltaT` (and `adt` but 
this is variable and not known in the recorder subroutine), then `tol` should be set relative 
to one of these parameters.


 1. Defining `tol = timeStamp * DBL_EPSILON / k`, where `k` is the allowable minimum ratio 
 of `adt/deltaT`, e.g. `k=1e-5`.
  - First precision error (late step) limit: `adt > deltaT * k`, error size: `adt`, 
  if `adt=deltaT` then 100% error in step. (working shown in section below)
  - Second precision error (early step) limit: `adt > DBL_EPSILON * timeStamp / k`, error size: `k * deltaT`

    Note that provided `k > 1e-7`, then second precision error does not occur because of the limit due to 
    the third precision error.

 2. Defining `tol = k * deltaT`
  - First precision error (late step) limit: `adt > DBL_EPSILON * timeStamp / k`, 
  error size: `adt`, if `adt=deltaT` then 100% error in step.
  - Second precision error (early step) limit: `adt > deltaT * k`, error size: `k * deltaT`

    Note that provided `k > 1e-7`, then first precision error does not occur because of the limit due 
    to the third precision error.

## Selected option

Option 2 is preferred since the error in the recorder can be limited to `k * deltaT`, 
provided that `k > 1e-7` and the third precision error is minimised. Ideally k would be 
close to `1e-7`, however, users may relax the third precision error limit, 
and therefore `1e-5` is appropriate, and still minimising the error size from the second 
precision error to reasonable levels comparable to the always persistent third precision error.


## Additional Documentation
$deltaT: time interval for recording. will record when next step is $deltaT greater than last 
recorder step. (optional, default: records at every time step). The reliable range for $deltaT is 
when the analysis time is greater than 1e-10 multiplied by the current time and greater than 1e-5 
multiplied $deltaT, assuming double precision tolerance of 1e-15.



## Working for Option 1 first precision error limit

`error < tol`

`deltaT/adt * DBL_EPSILON * timeStamp < timeStamp * DBL_EPSILON / k`

`adt > deltaT * k`


## Code examples

Example of precision error with no tolerance (current OpenSees release)

```python
# Implementation of recorder time limit in released version of OpenSees
import numpy as np
adt = 0.001
timeStamp = 0
nextTimeStampToRecord = 0
deltaT = 0.01
recorded_times = []
for i in range(1100):
    if timeStamp >= nextTimeStampToRecord:  # does not account for floating point precision
        recorded_times.append(timeStamp)
        nextTimeStampToRecord = timeStamp + deltaT
    timeStamp += adt
dts = np.diff(recorded_times)
print(np.sort(dts)[::-1][0])
# output = 0.0109999
assert max(dts) < deltaT + adt * 0.5, max(dts)  # This fails on MacOS 10.15
```

Example of third precision error

```python
# Implementation of analysis time step limit in released version of OpenSees
adt = 0.002
timeStamp = 1e12
initial_time = timeStamp
nextTimeStampToRecord = 0
deltaT = 0.01
recorded_times = []
n_analysis = 11000
for i in range(n_analysis):
    timeStamp += adt
expect_time_passed = adt * n_analysis
analysis_time_passed = timeStamp - initial_time
print(expect_time_passed, analysis_time_passed)
# 1.1, 2.1484375  # for MacOS 10.15
```

Example of previous commit failing using factor of 1e-12

```python
# Implementation of recorder time limit in previous commit of pull request
import numpy as np
adt = 0.0002
timeStamp = 10e4
nextTimeStampToRecord = 0
deltaT = 0.001
recorded_times = []
for i in range(110):
    # does not sufficiently account for cumulative error in timeStamp between recordings
    if timeStamp - nextTimeStampToRecord > -deltaT * 1e-12:
        recorded_times.append(timeStamp)
        nextTimeStampToRecord = timeStamp + deltaT
    timeStamp += adt
dts = np.diff(recorded_times)
print(np.sort(dts)[::-1][0])
# output 0.0011999999696854502
assert max(dts) < deltaT + adt * 0.5, max(dts)  # this fails
```

Example of option 1

```python
# Implementation of recorder time limit using option 1
# Analysis time step (adt) is set very close to upper limit for both
# floating point precision errors for MacOS 10.15 but still passes
# decreasing adt, increasing deltaT or increasing initial timeStamp will result
# in inaccurate recording.
import numpy as np
print(np.finfo(float).eps)
# output = 2.220446049250313e-16
eps = np.finfo(float).eps
adt = 0.000005
timeStamp = 2e5
nextTimeStampToRecord = 0
deltaT = 0.1
allowable_steps = 1e5
early_step_lim = timeStamp * eps * allowable_steps
late_step_lim = deltaT / allowable_steps
print('satisfies min_adt_lim_for_early_step: ', early_step_lim, adt, adt > early_step_lim)
print('satisfies min_adt_lim_for_late_step: ', late_step_lim, adt, adt > late_step_lim)
recorded_times = []
n_records = 4
n_steps = int(n_records * deltaT / adt)
for i in range(n_steps):
    equiv_zero_lim = -timeStamp * eps * allowable_steps
    if timeStamp - nextTimeStampToRecord > equiv_zero_lim:
        recorded_times.append(timeStamp)
        nextTimeStampToRecord = timeStamp + deltaT
    timeStamp += adt
dts = np.diff(recorded_times)
print('max_dt: ', max(dts))
print('min_dt: ', min(dts))
print('error: ', max(dts) - deltaT, min(dts) - deltaT)
print(np.sort(dts)[::-1][0])
# output = 0.10000017937272787
```

Example of option 2
```python
# Implementation of recorder time limit in latest commit of pull request
import numpy as np
adt = 0.0002
timeStamp = 10e8
nextTimeStampToRecord = 0
deltaT = 0.001
recorded_times = []
for i in range(110):
    if timeStamp - nextTimeStampToRecord > -deltaT * 1e-5:
        recorded_times.append(timeStamp)
        nextTimeStampToRecord = timeStamp + deltaT
    timeStamp += adt
dts = np.diff(recorded_times)
print(np.sort(dts)[::-1][0])
# output = 0.0010001659393310547
assert max(dts) < deltaT + adt * 0.5, max(dts)  # this passes
```
