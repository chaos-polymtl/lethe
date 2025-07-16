=============
Pull Requests
=============

Title
-----

Each pull request should only either fix an issue, add a single feature or restructure the code for one aspect, plus a test, if needed, to certify it.
If there are many aspects to change, it is better to create multiple branches and create multiple pull requests. The title should start with an action verb and only contain one idea.

Format
------

Lethe uses clang-format to have a uniform indentation across all source files. 
You should always run the indentation script in the contrib folder before creating a PR: ``contrib/utilities/indent-all``. 
If there are compatibility issues, you can run before that script either ``./contrib/utilities/download_clang_format`` or  ``./contrib/utilities/compile_clang_format``. 
For parameter files, there is a specific script that can be run as follows:  ``prmindent -i name_of_file.prm``.
Alternatively, you can run the script ``contrib/utilities/pre-commit.sh`` to indent both source and parameter files.

Review Process
--------------

* At least 2 reviewers for minor changes, like a bug fix or a small feature
* At least 3 reviewers for bigger changes like bigger features and architectural reconfiguration.
* The reviewer should be notified before opening the pull request. The reviewers should be selected from already existing contributor to the code. (See the README file).

Reviewers Responsibility
------------------------

* Give a review of the code implementation and general functionality of the code.
* Give a review on the code description and comments.
* Verify that the pull request meets the requirement described above regarding testing and the format of the pull request.
* Give the review in a timely manner.

Description
-----------

Without an explicit request to the person responsible for merging the PR, all commits will be squashed when merging. The pull request description will be the only remaining description of the changes left in the code history.
The description must therefore be filled in correctly using the template.
Depending on the type of changes, not all sections need to be filled in, but the contributor is expected to fill in the description as completely as possible.
This helps reviewers to understand what has been changed and facilitates future debugging.