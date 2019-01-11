import Pipeline.translations.wdl.wdl_parser as wdl

wdl_code = """
task my_task {
  File file
  String test = "default_value"
  command {
    ./my_binary --input=${file} > results
  }
  output {
    File result = "result"
    File results = "${test}"
  }
}

workflow my_wf {
  call my_task
}
"""


ww = wdl.parse(wdl_code).ast()

print(ww.dumps(indent=2))
body = ww.attr('body')
task = body[0]
print(task.attr('name').source_string)
task_command = task.attr('sections')[0]
print(task_command.attr("parts"))




