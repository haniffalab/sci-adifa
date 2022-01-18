# Contributing

Adifa is an open source project, and contributions are very welcome!

## Ways to contribute

Whether you're a developer, a data scientist or an academic, there are lots of ways to contribute. Here's a few ideas:

- [Install Adifa on your computer](https://haniffalab.com/sci-adifa/installation.html) and kick the tires. Does it work? Does it do what you'd expect? If not, [open an issue](https://github.com/haniffalab/sci-adifa/issues/new/choose) and let us know.
- Comment on some of the project's [open issues](https://github.com/haniffalab/sci-adifa/issues). Have you experienced the same problem? Know a work around? Do you have a suggestion for how the feature could be better?
- Read through the [documentation](https://haniffalab.com/sci-adifa/index.html), and click the "improve this page" button, any time you see something confusing, or have a suggestion for something that could be improved.
- Find an [open issue](https://github.com/haniffalab/sci-adifa/issues) (especially [those labeled `help-wanted`](https://github.com/haniffalab/sci-adifa/issues?q=is%3Aopen+is%3Aissue+label%3A%22help+wanted%22)), and submit a proposed fix. If it's your first pull request, we promise we won't bite, and are glad to answer any questions.
- Help evaluate [open pull requests](https://github.com/haniffalab/sci-adifa/pulls), by testing the changes locally and reviewing what's proposed.

## Code of Conduct

[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](code_of_conduct.md)

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation. This Code of Conduct is adapted from the [Contributor Covenant](http://contributor-covenant.org), version 2.1, available at [https://www.contributor-covenant.org/version/2/1/code_of_conduct/](https://www.contributor-covenant.org/version/2/1/code_of_conduct/)

## Before filing an issue

- Search the repository (also google) to see if someone has already reported the same issue. This allows contributors to spend less time responding to issues, and more time adding new features!
- Please provide a minimal complete verifiable example for any bug. If you're not sure what this means, check out [this blog post](http://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports) by Matthew Rocklin or [this definition](https://stackoverflow.com/help/mcve) from StackOverflow.
- Let us know about your environment.

## Submitting a pull request

### Pull requests generally

- The smaller the proposed change, the better. If you'd like to propose two unrelated changes, submit two pull requests.
- The more information, the better. Make judicious use of the pull request body. Describe what changes were made, why you made them, and what impact they will have for users.

### Submitting a pull request via github.com

Many small changes can be made entirely through the github.com web interface.

1. Navigate to the file within [`haniffalab/sci-adifa`](https://github.com/haniffalab/sci-adifa) that you'd like to edit.
2. Click the pencil icon in the top right corner to edit the file
3. Make your proposed changes
4. Click "Propose file change"
5. Click "Create pull request"
6. Add a descriptive title and detailed description for your proposed change. The more information the better.
7. Click "Create pull request"

That's it! You'll be automatically subscribed to receive updates as others review your proposed change and provide feedback.

### Submitting a pull request via Git command line

1. Fork the project by clicking "Fork" in the top right corner of [`haniffalab/sci-adifa`](https://github.com/haniffalab/sci-adifa).
2. Clone the repository locally `git clone git@github.com:haniffalab/sci-adifa.git`.
3. Create a new, descriptively named branch to contain your change ( `git checkout -b my-awesome-feature` ).
4. Hack away, add tests. Not necessarily in that order.
5. Make sure everything still passes by running `pytest` (see the [testing section](https://haniffalab.com/sci-adifa/installation.html#test) in the documentation)
6. Push the branch up ( `git push origin my-awesome-feature` ).
7. Create a pull request by visiting `https://github.com/<your-username>/sci-adifa` and following the instructions at the top of the screen.

## Proposing updates to the documentation

We want the Adifa documentation to be as helpful as possible. We've open-sourced our docs and we welcome any pull requests if you find it lacking.

### How to submit changes

You can find the Sphinx source fiels for our documentation in the [sphinx](https://github.com/haniffalab/sci-adifa/tree/main/sphinx) directory. See the section above, [submitting a pull request](#submitting-a-pull-request) for information on how to propose a change.

One gotcha, all pull requests should be directed at the `main` branch (the default branch).

## A thank you

Thanks! Working with Adifa should be fun. If you find any of this hard to figure out, let us know so we can improve our process or documentation!
