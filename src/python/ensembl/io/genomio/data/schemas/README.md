JSON Schema project: https://json-schema.org


To validate json using https://github.com/Julian/jsonschema:

Create pyenv:
```
export PYENV_ROOT=$HOME/pyenv
export PATH=$PYENV_ROOT/bin:$PATH

# pyenv install 3.7.4 # uncomment for the first time
export PYENV_VERSION=3.7.4

eval "$(pyenv init -)"
```

Install jsonschema:
```
pip install jsonschema
```

To validate JSON with schema just load it with python:
```
python -c 'import json; import sys; json.load(open(sys.argv[1]))' seq_region.json
``` 

Copy example from example section to example.json
```
python -c 'import sys; import json; \
    print("\n".join([ json.dumps(e) for e in json.load(open(sys.argv[1]))["examples"]]))' \
  seq_region.json > examples.json
```

And run validation:
```
cat examples.json |
  sed -n '2p' |
  jsonschema -i /dev/stdin seq_region.json && echo OK || echo FAIL
```
