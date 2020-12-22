# Janis CLI

> This document is still a work-in-progress. Please refer to `janis <command> -h` for the most comprehensive documentation.

The Janis CLI is an easy way to interact with Janis on the command line. 
This is installed automatically when you install `janis-pipelines` from git.

```bash
$ janis -h
```

You can see there are a lot of sub-commands:

- [`run`](#run): Run a Janis workflow
- [`prepare`](#prepare): Prepare a Janis workflow
- [`init`](#init): Initialise a Janis configuration
- [`translate`](#translate): Translate a janis workflow to CWL or WDL
- [`inputs`](#inputs): Generate an input job file for a janis workflow
- [`watch`](#watch): Watch an existing Janis workflow
- [`pause`](#pause): DEV: If something goes wrong, gracefully shut down Janis without ABORTing the workflow
- [`abort`](#abort): Abort a running Janis Workflow
- [`rm`](#rm): Remove the output directory and metadata
- [`metadata`](#metadata): Print all known metadata about a workflow
- [`spider`](#spider): Get information about a tool
- [`query`](#query): Search known workflows by some criteria
- [`cleanup`](#cleanup): Cleanup the central db of workflow runs
- [`rawquery`](#rawquery): Perform a raw SQL query on the sqlite database of a task
- [`wait`](#wait): Wait for all workflows to finish before returning, exiting with rc=3 if any of the workflows have failed.
- [`version`](#version): Print the versions of Janis and exit
- [`docs`](#docs): Attempts to open Janis documentation using `webbrowser`

## Run

Information about the `run` method

## Prepare

Information about the `prepare` method

## Init

Information about the `init` method

## Translate

Information about the `translate` method

## Inputs

Information about the `inputs` method

## Watch

Information about the `watch` method

## Resume

Information about the `resume` method

## Pause

Information about the `pause` method

## Abort

Information about the `abort` method

## Rm

Information about the `rm` method

## Metadata

Information about the `metadata` method

## Spider

Information about the `spider` method

## Query

Information about the `query` method

## Cleanup

Information about the `cleanup` method

## Rawquery

Information about the `rawquery` method

## Wait

Wait for all workflows to finish before returning

## Version

Prints all the versions of Janis components, similar to `janis -v`.
For example:

```bash
$ janis -v
# --------------------  -------
# janis-core            v0.11.0
# janis-assistant       v0.11.0
# janis-pipelines       v0.10.2
# janis-templates       v0.10.4
# janis-unix            v0.10.4
# janis-bioinformatics  v0.11.0
# --------------------  -------
```

## Docs

Attempts to open Janis documentation using a listed web browser, or an internal unix `webbrowser`.


## Other options

### Logging Levels

> The default log level is `INFO`

Janis has different levels of logs, and will only write to stderr the logs that are of "higher priority" that the current log level.

The log level params must be specified after `janis` but BEFORE the `command`, eg:

```bash
janis [LOG PARAMS] <method> [...options]
```

You can specify different log levels in the following ways:

- VERBOSE: `--verbose / --logVerbose`
- DEBUG: `-d / --debug / --logDebug`
- INFO: `--logInfo`
- WARN: `--logWarn`
- CRITICAL: `--logCritical`
- NONE (no logging): `--logNone`
- OR: `-L {VERB,DEBUG,INFO,WARN,CRITICAL,NONE}, --logLevel {VERB,DEBUG,INFO,WARN,CRITICAL,NONE}`
