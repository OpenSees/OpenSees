## Description
<!-- What does this PR do, and why? Link any related issue (e.g. "Fixes #1234"). -->

## Type of change
- [ ] Bug fix
- [ ] New element / material / section / other class
- [ ] Enhancement to existing code
- [ ] Build system / CI
- [ ] Documentation / examples / tests only

## AI code review (required)
Before opening this PR, please review your changes with an AI code-review tool
(for example Claude Code `/code-review`, GitHub Copilot code review, or similar).
See [CONTRIBUTING.md](https://github.com/OpenSees/OpenSees/blob/master/CONTRIBUTING.md#ai-code-review-before-you-open-a-pr) for how. We suggest starting with simply asking "Does this code contain errors?"

- [ ] I reviewed my changes with an AI code-review tool and fixed or answered the issues it raised

**Tool used:** <!-- e.g. Claude Code /code-review -->

**What it found and what I changed:**
<!-- A short summary. If it flagged something you chose NOT to change, say why. -->


## Checklist
- [ ] The PR contains only the changes it describes (no unrelated files, editor backups, build output, or local config)
- [ ] The code compiles with CMake on my platform: <!-- Windows / Linux / macOS -->
- [ ] Existing examples/tests that touch this code still give the same results (or the differences are explained above)

**If adding a new class (element, material, section, ...):**
- [ ] I have added a BSD-2 or BSD-3 Copyright notice to the code
- [ ] A test or example script is included (e.g. in `tests/` or `EXAMPLES/`)
- [ ] The class has a unique class tag in `SRC/classTags.h`
- [ ] The command is registered for both Tcl and Python interpreters
- [ ] `getCopy()`, `commitState()`, `revertToLastCommit()` and `revertToStart()` are implemented
- [ ] `sendSelf()` / `recvSelf()` are implemented and the class is added to `FEM_ObjectBrokerAllClasses.cpp` (needed for parallel runs)
- [ ] All new source files are added to the relevant `CMakeLists.txt`
- [ ] Documentation PR opened in [OpenSeesDocumentation](https://github.com/OpenSees/OpenSeesDocumentation) (link: )
