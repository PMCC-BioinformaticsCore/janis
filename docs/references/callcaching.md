# Call Caching

Call caching, or the ability for pipelines to use previous computed results. Each engine implements this separately, so this guide will assist you in setting up Janis, and how to configure the engines caching rules.

This feature can also be viewed as "resuming a workflow where you left off".

You may need to provide additional parameters when you run your workflow to ensure call-caching works. Please read this guide carefully.


## Configuring Janis

You need to include the following line in your [Janis Configuration](https://janis.readthedocs.io/en/latest/references/configuration.html#other-janis-options):

```yaml
# ...other options
call_caching_enabled: true
```

Strongly recommended:

- Use the `--development` run option as it ensures Cromwell will use the required database, and additionally sets the `--keep-intermediate-files` flag:

    - Remember, Janis will remove your execution directory (default: `<outputdir>/janis/execution`). And call caching only works if your intermediate files are immediately available. 

### CWLTool

No extra configuration should be required.

### Cromwell

> More information about how Cromwell call-caching works below

You're running on a local or shared filesystem (including HPCs), we strongly recommend running Cromwell version 50 or higher, and to [use the `fingerprint` hashing strategy ([Docs](https://cromwell.readthedocs.io/en/stable/Configuring/#local-filesystem-options) | [PR](https://github.com/broadinstitute/cromwell/pull/5450)):

```yaml
cromwell:
  call_caching_method: fingerprint
```

You MUST additionally run Janis with the `--mysql` (recommended: `--development`) flag, Cromwell [relies on a database](https://cromwell.readthedocs.io/en/stable/tutorials/PersistentServer/) for call-caching to work (unless you're running your own Cromwell server).

## How does call-caching work

### CWlTool

_Needs further clarification_.

### Cromwell


> More information: [Cromwell: Call Caching](https://cromwell.readthedocs.io/en/stable/cromwell_features/CallCaching/)

More docs are in progress (here is [our investigation](https://github.com/broadinstitute/cromwell/issues/5346)), but specifically, Cromwell hashes all the components of your task:

- output count
- runtime attributes
- output expression
- input count
- backend name
- command template
- input - hash of the file, dependent on the filesystem (local, GCS, S3) - see more below.

and determines a total hash for the call you are going to make. Cromwell will then check in its database to see if it's made a call that exactly matches this hash, and if so uses the stored result.

If this _call_ has NOT been cached, it will compute the result and then store it against the original hash for future use.


#### Input hash by filesystem

If you use a blob storage, like:

- Google Cloud Storage (GCS: `gs://`)
- Amazon Simple Storage Service (S3: `s3://`)

Cromwell can use the object id (the `etag`), as any modification to these files changes the object identifier.

The HTTP 

##### Local filesystem

In a local filesystem, you don't have this object ID luxury, so there are a number of hashing strategies:

- `file`: md5 hash of the file content (bad for large files),
- `path`: computes an md5 hash of the file path (doesn't work for containers),
- `path+modtime`: computes md5 hash of file path + modtime (doesn't work for containers),
- `xxhash` ([PROPOSED](https://github.com/broadinstitute/cromwell/pull/5450), > 50): Quicker hashing strategy for entire file. 
- `fingerprint` ([PROPOSED](https://github.com/broadinstitute/cromwell/pull/5450), > 50): Uses the `xxhash` strategy on (the first 10MB of file + mod time).

