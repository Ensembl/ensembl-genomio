## [Gitlab CI/CD](https://docs.gitlab.com/ee/ci/) pipelines
### Location [cicd/gitlab](https://github.com/Ensembl/ensembl-genomio/blob/main/cicd/gitlab)

Some Gitlab based CI/CD helper pipelines

### Description

### Setup
0) Import your GitHub project from GitLab.
See this [Gitlab documentation](https://docs.gitlab.com/ee/user/project/import/github.html).

1) Enable "Require status checks to pass before merging"? in "Settings"/"Branch protection rules" section for the "main" branch.

2) Get repo

3) Create (if you need)  cicd/gitlab folders, add configs.
We use [cicd/gitlab/dot.gitlab-ci.yml](https://github.com/Ensembl/ensembl-genomio/blob/main/cicd/gitlab/dot.gitlab-ci.yml) for this repo instead the default one (see below).

4) Edit settings for CI/CD
[General pipelines] expand
Set "CI/CD configuration file" to [cicd/gitlab/dot.gitlab-ci.yml](https://github.com/Ensembl/ensembl-genomio/blob/main/cicd/gitlab/dot.gitlab-ci.yml)
Click [Save changes] at the end of this small section

5) Enable runners
[General pipelines] expand
[Runners] expand
Pick "Shared runners"

6) Customize email settings
Like it's stated on the [official documentation](https://docs.gitlab.com/ee/user/project/integrations/pipeline_status_emails.html):

* Go to *Settings > Integrations > Pipeline status emails*.
* Edit *Recipients* (a comma-separated list of email addresses)
* Select *Notify only broken pipelines*
* Select the branches
* *Save*

## A few notes on style

* We suggest separating logic for running various parts into separate pipelines and using different `trigger` jobs to invoke these pipelines.
As for now we have [cicd/gitlab/parts/](https://github.com/Ensembl/ensembl-genomio/blob/main/cicd/gitlab/parts/) folder for these needs.

* Feel free to use external pipelines (as `trigger` jobs) or other parts from [GitLab templates](https://gitlab.com/gitlab-org/gitlab/-/tree/master/lib/gitlab/ci/templates). Though for pipelines, please, add 
```
  trigger:
    ...
    forward:
      pipeline_variables: false
      yaml_variables: false
...
```
to your trigger jobs whenever possible.

* We're trying to use `:` decided "namespace" qualifiers for stage and job names.

* Feel free to use template/generic jobs with names starting from `.` and `extends` keyword (see [`extends` description](https://docs.gitlab.com/ee/ci/yaml/#extends)


