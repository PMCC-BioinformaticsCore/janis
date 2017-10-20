import argparse

from wehiwdl import *;


def main( opts ):
    print("WeHi Pipeline Definition Demo!");
    #fac = cwltool.factory.Factory()

    #echo = fac.make("echo.cwl")
    #result = echo(message="foo")
    #wehiwdl.identify();


    wf = Workflow( label = "My Workflow" );
    wfd = WorkflowDocumentation( desc = "Workflow Demonstrator", author = "Mohammad Bhuyan");
    wf.setDocumentation( wfd  );

    ###########################  Example: Very Simmple ################################

    task = Task();
    tskd = TaskDocumentation( desc = "Copy task", author = "Mohammad Bhuyan");
    task.setDocumentation( tskd );

    ins = TaskInputSet();
    inp1 = TaskInputSpec( id = "in1", meta = { "type" : "FILE", "uri" : "http://www.s3.com/files/1.txt" } );
    ins.addInputSpec( inp1 )
    task.setInputSet( ins );

    outs = TaskOutputSet();
    out1 = TaskOutputSpec( id = "out1", meta = { "type" : "FILE", "uri" : "file://~/temp/1.txt" } );
    outs.addOutputSpec( out1 );
    task.setOutputSet( outs );


    #Will use templating engine. See Jinja (http://jinja.pocoo.org/)
    cmd = Command( template = "cp <inputs['in1']> <outputs['out']>");
    task.setCommand( cmd );

    wf.addTask( step="COPY", task= task );

    # 1. Will use concept of registry of translator / emitter. "CWL" translator regisered.
    # 2. Will validate against schema and raise error
    wf.generateDefinition( format = "CWL", meta = { "file" : "~/temp/wf.cwl" } );


    ###########################  Example: Input Output Mapping ################################


    wfinputs = WorkflowInputSet();
    wfinput = WorkflowInputSpec( id = "win1", meta = { "type" : "FILE", "uri" : "http://www.s3.com/files/1.txt" } );
    wfinputs.addInputSpec( wfinput );
    wf.setInputSet( wfinputs );

    wfoutputs = WorkflowOutputSet();
    wfoutput = WorkflowOutputSpec( id = "wout1", meta = { "type" : "FILE", "uri" : "file://~/temp/1.txt" } );
    wfoutputs.addOutputSpec( wfoutput );
    wf.setOutputSet( wfinputs );

    task = Task();
    tskd = TaskDocumentation( desc = "Copy task", author = "Mohammad Bhuyan");
    task.setDocumentation( tskd );

    ins = TaskInputSet();
    inp1 = TaskInputSpec( id = "in1", meta = { "ref" : "win1" } ); # <--- Ref to workflow input
    ins.addInputSpec( inp1 )
    task.setInputSet( ins );

    outs = TaskOutputSet();
    out1 = TaskOutputSpec( id = "out1", meta = { "ref" : "wout1" } ); # <--- Ref to workflow output
    outs.addOutputSpec( out1 );
    task.setOutputSet( outs );

    #Will use templating engine. See Jinja (http://jinja.pocoo.org/)
    cmd = Command( template = "cp <inputs['in1']> <outputs['out']>");
    task.setCommand( cmd );

    wf.addTask( step="COPY", task= task );

    # 1. Will use concept of registry of translator / emitter. "CWL" translator regisered.
    # 2. Will validate against schema and raise error
    wf.generateDefinition( format = "CWL", meta = { "file" : "~/temp/wf.cwl" } );


    ###########################  Example: Nested Workflow ################################








if __name__ == "__main__":
    argprsr = argparse.ArgumentParser();
    #argprsr.add_argument("--config", required=False);
    opts = argprsr.parse_args()
    main( opts );
