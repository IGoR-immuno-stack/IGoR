## test suite description


## mechanism

- on github.com, create labels on the form `test::XXX` where XXX is the name of the test suite to runner
  (eg: test::unit)

- associated one or several os these labels to a PR/MR

- create pixi tasks names `test_XXX` which correspond to these tests

- the `resolve_test_selection.py` script is run by the CI actions and will detects the associated labels to each PR/MR and run the test suite accordingly


## local debug

### Create a JSON file to simulate `GITHUB_EVENT_PATH`

cat > /tmp/test-event.json << EOF
{
  "pull_request": {
    "labels": [
      {"name": "test::unit"},
      {"name": "test::integration"},
      {"name": "WIP"}
    ]
  }
}
EOF

### debug the previous simulated behavior

```
GITHUB_EVENT_PATH=/tmp/test-event.json pixi run resolve-test-selection --dry-run --default-task test
```

### run the tests effectively

```
GITHUB_EVENT_PATH=/tmp/test-event.json pixi run resolve-test-selection --execute --default-task test
```
