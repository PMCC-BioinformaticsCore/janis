import argparse
import atexit

# Import library of already defined inputs, steps and outputs

from pipeline_definition.pipeline_translator import PipelineTranslator

# Extend the translation system with user defined


def main(options):

  def at_exit():
    print("BYE!!")

  atexit.register(at_exit)

  # Get specified file to translate
  pdfile = options.pdfile

  # pdfile = "pd_1.yml"

  pd_translator = PipelineTranslator()
  pd_translator.translate(pdfile, outfile="/Users/mohammadbhuyan/Temp/out.pdx", overwrite_outfile=True)
  # pdx.translate(pdfile, outfile="/Users/mohammadbhuyan/Temp/out.pdx")


if __name__ == "__main__":
  argprsr = argparse.ArgumentParser()
  argprsr.add_argument('pdfile', help='Pipeline Definition file.')
  opts = argprsr.parse_args()
  main(opts)
