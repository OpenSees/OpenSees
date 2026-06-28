# OpenSees Source Code Repository

This git repository contains all revisions to OpenSees source code since Version 2.3.2.

Older revisions of the code are available upon request.

If you plan on collaborating or even using OpenSees as your base code it is highly recommended that
you FORK this repo to your own account and work on it there. We will not allow anybody to write to
this repo. Only PULL requests will be considered. To fork the repo click on the FORK at the top of this page.

For a brief outline on forking we suggest:
https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow

For a brief introduction to using your new FORK we suggest:
https://www.atlassian.com/git/tutorials/saving-changes

## Documentation
The documentation for OpenSees is being moved to a parellel github repo:
https://github.com/OpenSees/OpenSeesDocumentation

The documentation (in its present form) can be viewed in the browser using the following url:
https://OpenSees.github.io/OpenSeesDocumentation

## Build Instructions
Steps to build OpenSees on Windows, Linux, and Mac:
https://opensees.github.io/OpenSeesDocumentation/developer/build.html

### Optional Sanitizer Builds

Developers using GNU or Clang compilers can enable runtime instrumentation while configuring CMake:

```console
cmake -S . -B build/asan -DOPS_ENABLE_ASAN=ON
cmake -S . -B build/ubsan -DOPS_ENABLE_UBSAN=ON
```

These options are disabled by default and are intended for debugging memory and undefined-behavior issues.

## Modeling Questions
Issues related to modeling questions will be closed. Instead, post your modeling questions on the OpenSees 
message board or in the OpenSees Facebook group.
+ https://opensees.berkeley.edu/community
+ https://facebook.com/groups/opensees
