## [Gitlab CI/CD](https://docs.gitlab.com/ee/ci/) pipelines
### Location [gitlab/cicd](../gitlab/cicd)

Some Gitlab based CI/CD helper pipelines

### Description

### Setup
0) Import your GitHub project from GitLab.
See this [Gitlab documentation](https://docs.gitlab.com/ee/user/project/import/github.html).

1) Enable "Require status checks to pass before merging"? in "Settings"/"Branch protection rules" section for the "main" branch.

2) Get repo

3) Create (if you need)  cicd/gitlab folders, add configs.
We use [cicd/gitlab/dot.gitlab-ci.yml](../cicd/gitlab/dot.gitlab-ci.yml) for this repo instead the default one (see below).

4) Edit settings for CI/CD
[General pipelines] expand
set "CI/CD configuration file" to [cicd/gitlab/dot.gitlab-ci.yml](../cicd/gitlab/dot.gitlab-ci.yml)
click [Save changes] at the end of this small section

5) enable runners
General pipelines] expand
[Runners] expand
pick "Shared runners"

## A few notes on style

* We suggest separating logic for running various parts into separate pipelines and using different `trigger` jobs to invoke these pipelines.
As for now we have [cicd/gitlab/parts/](../cicd/gitlab/parts/) folder for these needs.

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


