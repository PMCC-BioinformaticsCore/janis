import Pipeline.translations.wdl.wdl_parser as wdl

wdl_code = """
task my_task {
  File file
  command {
    ./my_binary --input=${file} > results
  }
  output {
    File results = "results"
  }
}

workflow my_wf {
  call my_task
}
"""


ww = wdl.parse(wdl_code).ast()

print(ww.dumps(indent=2))
print(ww.attr('body')[1].attr('name').source_string)



