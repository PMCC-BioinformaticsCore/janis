import Pipeline.translations.wdl.wdl_parser as wdl

wdl_code = """
task hello {
  String salutation
  String name

  command {
    echo '${salutation}, ${name}!'
  }
  output {
    String response = read_string(stdout())
  }
}

workflow test {
  String greeting
  call hello {
    input: salutation=greeting
  }
  call hello as hello2 {
    input: salutation=greeting + " and nice to meet you"
  }
}
"""




ww = wdl.parse(wdl_code).ast()

print(ww.dumps(indent=2))
body = ww.attr('body')
task = body[0]
print(task.attr('name').source_string)
task_command = task.attr('sections')[0]
print(task_command.attr("parts"))


