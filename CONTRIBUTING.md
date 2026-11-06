# Contributing to OpenSees

Thank you for contributing to OpenSees. This page explains how to get a pull
request (PR) ready so it can be reviewed and merged quickly.

> **Modeling questions are not handled here.** Please post them on the
> [OpenSees message board](https://opensees.berkeley.edu/community)

## Workflow

1. **Fork** this repository to your own GitHub account. Nobody pushes directly
   to this repo; only pull requests are considered.
2. Create a branch in your fork for your change.
3. Build and test your change (see below).
4. We sugget **Running an AI code review** on your changes (see below) and fix what it finds.
5. Open a pull request against `master` and fill in the PR template.

## AI code review before you open a PR

OpenSees maintainers are a small team reviewing a large number of PRs. We ask
every contributor to review their own changes with an AI code-review tool
before opening a PR. These tools are good at finding the kinds of bugs that
take reviewers the most time, such as:

- memory leaks, missing `delete`, and double frees
- uninitialized member variables
- trial and committed state not handled correctly in `commitState()`,
  `revertToLastCommit()` or `revertToStart()`
- `getCopy()` not copying all of the state
- `sendSelf()` / `recvSelf()` not sending or receiving the same data in the same order
- array index errors and divide-by-zero
- input parsing that does not check the number or type of arguments

You can use any tool you like. In the PR template, tell us which tool you used,
what it found, and what you changed. If you chose not to fix something it
flagged, say why. That is fine; the AI is not always right.

### Example: Claude Code

1. Install [Claude Code](https://docs.claude.com/en/docs/claude-code/overview).
2. In your local clone, check out your branch and start it:
   ```bash
   cd OpenSees
   git checkout my-feature-branch
   claude
   ```
3. Ask for a review of your branch's changes against `master`:
   ```
   /code-review
   ```
   or, in plain language:
   ```
   Review the changes on this branch compared to master for bugs. This is
   OpenSees; check state handling (commit/revert), getCopy, sendSelf/recvSelf,
   memory management, and input argument checking.
   ```
4. Fix the issues, re-run the review, and then open your PR.

### Other tools

GitHub Copilot code review, Cursor, ChatGPT/Codex, Gemini and others also work.
If your tool needs you to paste code, paste your diff (`git diff master...HEAD`)
along with the prompt in step 3 above.

**Note:** an AI review does not replace testing. You are still responsible for
building the code and checking that it gives correct results.

**Note:** we will probably use AI to test what it is submuitted, assuming it compiles. Code that fails to compile will **NOT** be pulled in. We will not fix the compilation issues. Once the PR is made, if you have submitted code before automatic testing will be performed and the PR page will indicate it it compiled. If new, before we accept the pull request,  we will run the tests in github and notify you if the compilation fails.

## What makes a good PR

- **Keep it focused.** One feature or fix per PR. Do not include unrelated
  formatting changes, editor backup files, build output, or local config files.
- **It must build.** Make sure it compiles with CMake. CI builds every PR.
- **Include a test.** For new code, add a test or example script that shows it
  working, e.g. a `pytest` script in [`tests/`](tests/README.md) or a Tcl/Python
  example in `EXAMPLES/`.
- **Do not change existing results** without explaining why in the PR.

## Adding a new element, material, section, ...

### Required:

Most new code in OpenSees is a new class. Please make sure it:
- has a BSD-2 or BSD-3 Copyright notice (you do not need to copy the existing OpenSees copyright notice .. ity just must be BSD-2 or BSD-3)
- has a unique class tag in `SRC/classTags.h`
- implements `getCopy()`, `commitState()`, `revertToLastCommit()` and `revertToStart()`
- implements `sendSelf()` / `recvSelf()`, and is added to
  `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` so it works in
  parallel runs
- has its source files added to the relevant `CMakeLists.txt`
- has documentation: open a PR in
  [OpenSeesDocumentation](https://github.com/OpenSees/OpenSeesDocumentation)

### Helpful:

- is registered for **both** the Tcl and Python interpreters
  (for uniaxial materials, see `SRC/interpreter/OpenSeesUniaxialMaterialCommands.cpp`
  and `SRC/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp`)

