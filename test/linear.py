import unittest

from pipeline_definition.pipeline_translator import PipelineTranslator
import json
import examples.bio_informatics


_yml = """
inputs:
  fastq:
    SequenceReadArchivePaired:
      forward-pattern: '*_R1.fastq.gz'
      backward-pattern: '*_R2.fastq.gz'
  ref:
    BAM:
      path: 'path/to/reference'

steps:
  - step1:
      trim:
        trimmer : 'trimmomatic'
  - step2:
      input_scope: ['ref']
      align:
        aligner: 'bwa'
  - step3:
      dedup:
  - step4:
      input_scope: ['step2', 'step3']
      call:
"""

_expected = json.loads("""
{
    "workflow": {
        "inputs": {
            "fastq": {
                "type": "SequenceReadArchivePaired"
            },
            "ref": {
                "type": "BAM"
            }
        },
        "flow": {
            "untagged": {
                "steps": {
                    "1": {
                        "step": "step1",
                        "type": "trim",
                        "step-inputs": {
                            "read": {
                                "type": "SequenceReadArchivePaired",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "fastq",
                                            "type": "SequenceReadArchivePaired",
                                            "step": "input-step",
                                            "tag": "input"
                                        }
                                    }
                                }
                            }
                        },
                        "step-outputs": {
                            "trimmed": {
                                "type": "SequenceReadArchivePaired"
                            }
                        }
                    },
                    "2": {
                        "step": "step2",
                        "type": "align",
                        "step-inputs": {
                            "read": {
                                "type": "SequenceReadArchivePaired",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "trimmed",
                                            "type": "SequenceReadArchivePaired",
                                            "step": "step1",
                                            "tag": "untagged"
                                        },
                                        "2": {
                                            "id": "fastq",
                                            "type": "SequenceReadArchivePaired",
                                            "step": "input-step",
                                            "tag": "input"
                                        }
                                    }
                                }
                            },
                            "reference": {
                                "type": "REFERENCE",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "ERROR": "Failed to find any candidate!!!!!"
                                    }
                                }
                            }
                        },
                        "step-outputs": {
                            "alignedbamfile": {
                                "type": "BAM"
                            }
                        }
                    },
                    "3": {
                        "step": "step3",
                        "type": "dedup",
                        "step-inputs": {
                            "bamfile": {
                                "type": "BAM",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "alignedbamfile",
                                            "type": "BAM",
                                            "step": "step2",
                                            "tag": "untagged"
                                        },
                                        "2": {
                                            "id": "ref",
                                            "type": "BAM",
                                            "step": "input-step",
                                            "tag": "input"
                                        }
                                    }
                                }
                            }
                        },
                        "step-outputs": {
                            "bamfile": {
                                "type": "BAM"
                            }
                        }
                    },
                    "4": {
                        "step": "step4",
                        "type": "call",
                        "step-inputs": {
                            "alignedbamfile": {
                                "type": "BAM",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "bamfile",
                                            "type": "BAM",
                                            "step": "step3",
                                            "tag": "untagged"
                                        },
                                        "2": {
                                            "id": "alignedbamfile",
                                            "type": "BAM",
                                            "step": "step2",
                                            "tag": "untagged"
                                        }
                                    }
                                }
                            },
                            "reference": {
                                "type": "REFERENCE",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "ERROR": "Failed to find any candidate!!!!!"
                                    }
                                }
                            }
                        },
                        "step-outputs": {
                            "read": {
                                "type": "VCF"
                            }
                        }
                    }
                }
            }
        }
    }
}""")


class LinearPipeline(unittest.TestCase):

  def test_graph(self):
    translator = PipelineTranslator(debug=True)
    translator.translate_string(_yml)
    tr_json = json.loads(translator.pipeline())
    self.assertTrue(tr_json == _expected)


if __name__ == '__main__':
    unittest.main()
