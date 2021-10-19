# Contribute to DIIM

DIIM provides the dynamic inoperability input-output model. This document 
describes how you can contribute.

## How can I contribute?

### Reporting bugs

You can report bugs by [open up a new issue](https://github.com/stigrs/diim/issues/new) 
on GitHub. Please ensure that the bug not already reported by searching under
[Issues](https://github.com/stigrs/diim/issues). Be sure to include a *title 
and clear description*, as much relevant information as possible, and a *code 
sample* or an *executable test case* demonstrating the expected behavior that is 
not occurring. 

### Fixing bugs

Open a new GitHub pull request with the patch. Ensure that the pull request
description clearly describes the problem and the solution. Include the relevant
issue number if applicable. Before submitting, please read the coding style 
guide to know more about coding conventions.

### Suggesting enhancements

Suggest your change and start writing code in accordance with the coding 
style guide if you get positive feedback from the project owner. Propose the
change by forking the repo and open a pull request:
1. Create a fork
2. Clone the fork
3. Create a feature branch
4. Commit to the feature branch
5. Push the feature branch to the fork
6. Open a pull request, making sure that:
    * Unfinished work is not submitted
    * The commits are atomic and the commit messages are useful
    * The new code is tested and validated
    * The new code is documented

## Coding Style Guide
* The source code should follow the [C++ Core Guidelines](http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines),
especially when it comes to naming rules
* The source code shall follow the style supplied in the .clang-format file
* The source code shall be documented
* The source code shall be tested and validated with [Catch2](https://https://github.com/catchorg/Catch2)
test cases
