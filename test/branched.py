import unittest

from pipeline_definition.pipeline_translator import PipelineTranslator
import json
import examples.bio_informatics

_yml = """
inputs:
    fastq:
        type: SequenceReadArchivePaired
        forward-pattern: '*_R1.fastq.gz'
        backward-pattern: '*_R2.fastq.gz'
    ref:
        type: REFERENCE:
        path: 'path/to/reference'

steps:
    qc:
      tool: fastqc      # outputs REFERENCE
      input: fastq
    remove-adaptors:
      tool: trim        # outputs TrimmedReads
      reads2: fastq 
    align-to-human:
      tool: align       # outputs: BAMFILE
      trimmed-reads: remove-adaptors/output
      reference: qc/output
    dedup:
        tool: dedup     # outputs BAMFILE
        bamfile: align-to-human/output
    intersect-genic:
        # input: [dedup]
        read: trim
        tool: bedtools-intersect
        split: true
    intersect-nongenic:
        # input: [dedup]
        read: trim
        tool: bedtools-intersect
        reportNoOverlaps: true
"""

_expected = json.loads("""
{
    "workflow": {
        "inputs": {
            "fastq": {
                "type": "SequenceReadArchivePaired"
            },
            "ref": {
                "type": "REFERENCE"
            }
        },
        "flow": {
            "untagged": {
                "steps": {
                    "1": {
                        "step": "qc",
                        "type": "fastqc",
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
                            "reports": {
                                "type": "Text"
                            }
                        }
                    },
                    "2": {
                        "step": "remove-adaptors",
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
                    "3": {
                        "step": "align-to-human",
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
                                            "step": "remove-adaptors",
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
                                        "1": {
                                            "id": "ref",
                                            "type": "REFERENCE",
                                            "step": "input-step",
                                            "tag": "input"
                                        }
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
                    "4": {
                        "step": "dedup",
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
                                            "step": "align-to-human",
                                            "tag": "untagged"
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
                    "5": {
                        "step": "intersect-genic",
                        "type": "bedtools-intersect",
                        "step-inputs": {
                            "read": {
                                "type": "SequenceReadArchivePaired",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "trimmed",
                                            "type": "SequenceReadArchivePaired",
                                            "step": "remove-adaptors",
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
                            }
                        },
                        "step-outputs": {
                            "reports": {
                                "type": "Text"
                            }
                        }
                    },
                    "6": {
                        "step": "intersect-nongenic",
                        "type": "bedtools-intersect",
                        "step-inputs": {
                            "read": {
                                "type": "SequenceReadArchivePaired",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "trimmed",
                                            "type": "SequenceReadArchivePaired",
                                            "step": "remove-adaptors",
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
                            }
                        },
                        "step-outputs": {
                            "reports": {
                                "type": "Text"
                            }
                        }
                    }
                }
            }
        }
    }
}""")


class BranchedPipeline(unittest.TestCase):

  def test_graph(self):
    translator = PipelineTranslator(debug=True)
    translation = translator.translate_string(_yml)
    tr_json = json.loads(translation)
    self.assertTrue(tr_json == _expected)


if __name__ == '__main__':
    unittest.main()


