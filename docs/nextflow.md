# Nextflow related documentation

## Installation

Here's as simple bash snippet that allows you to install nextflow into a dedicate folder without keeping any bits of the Nexflow (`.nextflow` folder) in your HOME dir.
We using bash here, though feel free to use any other shell your prefer.

If, for any reason, you're not happy with the suggested approach, feel free to use the default Nextflow installation instructions from [https://www.nextflow.io/index.html#GetStarted](https://www.nextflow.io/index.html#GetStarted).

We suggest following:

Install, tweaking that `HOME` env variable (`NEXTFLOW_HOME` env variable intoduced):

```
# add NEXTFLOW_HOME env
export NEXTFLOW_HOME=$(pwd) # or whatever
echo NEXTFLOW_HOME "$NEXTFLOW_HOME"

# get nextflow and install almost like here: https://www.nextflow.io/index.html#GetStarted
wget -O - https://get.nextflow.io  > nextflow.install.bash
# review and run, tweaking $HOME
cat nextflow.install.bash | HOME="${NEXTFLOW_HOME}" bash -i 2>&1 | tee nextflow.install.log

# make `.nextflow` a bit more visible
mv .nextflow dot.nextflow
ln -s dot.nextflow .nextflow
```

Now, hopefully, you have `nextflow` installed, but we need to tweak runnable to make it use
of our `NEXTFLOW_HOME` env var. We'll rename the original binary and create a shell wrapper
to run it, substituting `HOME` env variable with the `NEXTFLOW_HOME` value. 

N.B. This wrapper can try coping up-to-date `~/.bash_profile` and `~/.bashrc` from `$HOME` into
`$NEXTFLOW_HOME` if you need this, i.e. for LSF spawned jobs or whatever. For this to happen, uncomment corresponding lines.

```
# create a wrapper to deal with the default location of `.nextflow`
mv nextflow{,.bin}

cat > nextflow << 'EOF'
# copy .bashrc / .bash_profile here (for using on the LSF nodes i.e.)
#   N.B. uncomment the following 2 lines if it's secure and you need this
# ( cp ~/.bashrc "${NEXTFLOW_HOME}" )
# ( cp ~/.bash_profile "${NEXTFLOW_HOME}" )

# tweak HOME env
export HOME="${NEXTFLOW_HOME}"
"${NEXTFLOW_HOME}"/nextflow.bin "$@"
EOF

chmod a+rx "${NEXTFLOW_HOME}"/nextflow
```

Run test as suggestet by https://www.nextflow.io/index.html#GetStarted
```
# run test, see https://www.nextflow.io/index.html#GetStarted
./nextflow run hello
```

