# Adding a tool to a toolbox

> This section is optional.

This section requires that you have [Git](https://git-scm.com/) configured, and a [GitHub](https://help.github.com/en/github/getting-started-with-github/signing-up-for-a-new-github-account) account.

## Fork

> A fork is a copy of a repository. Forking a repository allows you to freely experiment with changes without affecting the original project. - GitHub: [Fork a Repo](https://help.github.com/en/github/getting-started-with-github/fork-a-repo)

> Further reading: [Fork an example repository](https://help.github.com/en/github/getting-started-with-github/fork-a-repo#fork-an-example-repository) 

If you want to add tools to an existing toolbox, you could for example, 

- Fork janis-unix.
- Add and test new tools,
- Create a [Pull Request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request#creating-the-pull-request) to integrate your changes.

To create a fork of `janis-unix`:

1. Navigate to the Janis-Unix repository:
    ```
    https://github.com/PMCC-BioinformaticsCore/janis-unix
    ```

2. In the top-right corner of the page, click **Fork**.

    ![Fork a repository](https://help.github.com/assets/images/help/repository/fork_button.jpg)

### Clone your fork

```bash
# 1. Clone your repoistory
git clone https://github.com/yourusername/janis-unix.git

# 2. Go into the checked out code
cd janis-bioinformatics

# 3. Add the original repo as 'upstream'
git remote add upstream https://github.com/PMCC-BioinformaticsCore/janis-unix.git
```

#### Keeping your fork up to date

> Further reading: [Syncing a fork](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/syncing-a-fork)

```bash
# 1. Get the changes from the original changes
git fetch upstream master

# 2. Merge these changes into your current branch
git merge upstream/master

# 3. Push these updates to your fork
git push
```

You may have to deal with conflicts if changes in the upstream have modified the same files you've locally modified.


## Making a change

First off, make sure your fork is up to date. Then let's checkout to a new branch:

```bash
git checkout -b add-greet-tool
```

Our new tool is going to simply print `"Hello, " + $name` to the console.

```bash
vim janis_unix/tools/greet.py
```

```python
from janis_core import String, ToolInput, ToolOutput, ToolArgument, Stdout, Boolean
from .unixtool import UnixTool


class Greet(UnixTool):
    def tool(self):
        return "greet"

    def friendly_name(self):
        return "Greet"

    def base_command(self):
        return "echo"

    def arguments(self):
        return [ToolArgument("Hello, ", position=0)]

    def inputs(self):
        return [ToolInput("name", String(), position=1)]

    def outputs(self):
        return [ToolOutput("out", Stdout())]
```

Let's make sure we add the Greet tool to be exported by this toolbox:

```bash
vim janis_unix/tools/__init__.py
```

```python
# ...other imports
from .hello import HelloWorkflow
from .greet import Greet
```


### Install the repository

The Janis environment at Peter Mac has been configured to allow your locally installed modules to override the modules from within the python environment.

You just have to reference the version of `pip3` outside of the Janis Python env: `/usr/bin/pip3`.

> Don't forget the `.` at the end of the pip install!

```bash
# --no-dependencies flag ensures we don't install janis core / assistant
/usr/bin/pip3 install --no-dependencies --user .
```

### Testing the change

Let's test the tool we just added:

```bash
janis spider Greet
```

_Voila!_ We added a tool to the janis-unix toolbox.

## Commit and push your changes

Let's commit our changes!

```bash
# 0. Check the status of your changes
git status

# 1. Stage our changes to commit (all of them)
git add .

# 2. Commit with a message
git commit -m "Adds new tool 'Greet'"
```

We want to push your branch `add-greet-tool` to the remote to create a pull request.

```
# Set the upstream on our branch, and push to remote 'origin'
git push --set-upstream origin add-greet-tool

Total 0 (delta 0), reused 0 (delta 0)
remote: 
remote: Create a pull request for 'add-greet-tool' on GitHub by visiting:
remote:      https://github.com/yourname/janis-unix/pull/new/add-greet-tool
remote: 
To github.com:yourname/janis-unix.git
 * [new branch]      add-greet-tool -> add-greet-tool
Branch 'add-greet-tool' set up to track remote branch 'add-greet-tool' from 'origin'.
```

Create a pull request by clicking the link:
- https://github.com/yourname/janis-unix/pull/new/add-greet-tool