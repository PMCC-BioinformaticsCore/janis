import unittest

from pipeline_definition.pipeline_translator import PipelineTranslator
import json
import examples.bio_informatics

_yml = """
inputs:
  tumour:
    SequenceReadArchivePaired:
      forward-pattern: '*_R1.fastq.gz'
      backward-pattern: '*_R2.fastq.gz'
  normal:
    SequenceReadArchivePaired:
      forward-pattern: '*_R1.fastq.gz'
      backward-pattern: '*_R2.fastq.gz'
  ref:
    reference:
      path: 'path/to/reference'

steps:
  - trim_tumour:
      tag: 'tumour'
      trim:
        trimmer : 'trimmomatic'
  - align_tumour:
      tag: 'tumour'
      align:
        aligner: 'bowtie2'
  - sort_tumour:
      tag: 'tumour'
      sort:
  - call_tumour:
      tag: 'tumour'
      call:
  - trim_normal:
      tag: 'normal'
      trim:
        trimmer : 'trimmomatic'
  - align_normal:
      tag: 'normal'
      align:
        aligner: 'bowtie2'
  - sort_normal:
      tag: 'normal'
      sort:
  - call_normal:
      tag: 'normal'
      call:
  - detect_mutation:
      joint_call:
        caller: mutect
        tumour_tag: '#tumour'
        normal_tag: '#normal'

"""

_expected = json.loads("""
{
    "workflow": {
        "inputs": {
            "tumour": {
                "type": "SequenceReadArchivePaired"
            },
            "normal": {
                "type": "SequenceReadArchivePaired"
            },
            "ref": {
                "type": "reference"
            }
        },
        "flow": {
            "tumour": {
                "steps": {
                    "1": {
                        "step": "trim_tumour",
                        "type": "trim",
                        "step-inputs": {
                            "read": {
                                "type": "SequenceReadArchivePaired",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "tumour",
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
                                "type": "TrimmedReads"
                            }
                        }
                    },
                    "2": {
                        "step": "align_tumour",
                        "type": "align",
                        "step-inputs": {
                            "trimmed reads": {
                                "type": "TrimmedReads",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "trimmed",
                                            "type": "TrimmedReads",
                                            "step": "trim_tumour",
                                            "tag": "tumour"
                                        }
                                    }
                                }
                            },
                            "reference": {
                                "type": "reference",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "ref",
                                            "type": "reference",
                                            "step": "input-step",
                                            "tag": "input"
                                        }
                                    }
                                }
                            }
                        },
                        "step-outputs": {
                            "alignedbamfile": {
                                "type": "bam"
                            }
                        }
                    },
                    "3": {
                        "step": "sort_tumour",
                        "type": "sort",
                        "step-inputs": {
                            "bamfile": {
                                "type": "bam",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "alignedbamfile",
                                            "type": "bam",
                                            "step": "align_tumour",
                                            "tag": "tumour"
                                        }
                                    }
                                }
                            }
                        },
                        "step-outputs": {
                            "sortedbamfile": {
                                "type": "sortedbam"
                            }
                        }
                    },
                    "4": {
                        "step": "index_tumour",
                        "type": "index",
                        "step-inputs": {
                            "sortedbamfile": {
                                "type": "bam",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "sortedbamfile",
                                            "type": "sortedbam",
                                            "step": "sort_tumour",
                                            "tag": "tumour"
                                        },
                                        "2": {
                                            "id": "alignedbamfile",
                                            "type": "bam",
                                            "step": "align_tumour",
                                            "tag": "tumour"
                                        }
                                    }
                                }
                            }
                        },
                        "step-outputs": {
                            "indexedfile": {
                                "type": "bamindex"
                            }
                        }
                    },
                    "5": {
                        "step": "call_tumour",
                        "type": "call",
                        "step-inputs": {
                            "indexedbam": {
                                "type": "bam",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "alignedbamfile",
                                            "type": "bam",
                                            "step": "align_tumour",
                                            "tag": "tumour"
                                        }
                                    }
                                }
                            },
                            "reference": {
                                "type": "reference",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "ref",
                                            "type": "reference",
                                            "step": "input-step",
                                            "tag": "input"
                                        }
                                    }
                                }
                            }
                        },
                        "step-outputs": {
                            "calls": {
                                "type": "vcf"
                            }
                        }
                    }
                }
            },
            "normal": {
                "steps": {
                    "1": {
                        "step": "trim_normal",
                        "type": "trim",
                        "step-inputs": {
                            "read": {
                                "type": "SequenceReadArchivePaired",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "normal",
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
                                "type": "TrimmedReads"
                            }
                        }
                    },
                    "2": {
                        "step": "align_normal",
                        "type": "align",
                        "step-inputs": {
                            "trimmed reads": {
                                "type": "TrimmedReads",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "trimmed",
                                            "type": "TrimmedReads",
                                            "step": "trim_normal",
                                            "tag": "normal"
                                        }
                                    }
                                }
                            },
                            "reference": {
                                "type": "reference",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "ref",
                                            "type": "reference",
                                            "step": "input-step",
                                            "tag": "input"
                                        }
                                    }
                                }
                            }
                        },
                        "step-outputs": {
                            "alignedbamfile": {
                                "type": "bam"
                            }
                        }
                    },
                    "3": {
                        "step": "sort_normal",
                        "type": "sort",
                        "step-inputs": {
                            "bamfile": {
                                "type": "bam",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "alignedbamfile",
                                            "type": "bam",
                                            "step": "align_normal",
                                            "tag": "normal"
                                        }
                                    }
                                }
                            }
                        },
                        "step-outputs": {
                            "sortedbamfile": {
                                "type": "sortedbam"
                            }
                        }
                    },
                    "4": {
                        "step": "index_normal",
                        "type": "index",
                        "step-inputs": {
                            "sortedbamfile": {
                                "type": "bam",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "sortedbamfile",
                                            "type": "sortedbam",
                                            "step": "sort_normal",
                                            "tag": "normal"
                                        },
                                        "2": {
                                            "id": "alignedbamfile",
                                            "type": "bam",
                                            "step": "align_normal",
                                            "tag": "normal"
                                        }
                                    }
                                }
                            }
                        },
                        "step-outputs": {
                            "indexedfile": {
                                "type": "bamindex"
                            }
                        }
                    },
                    "5": {
                        "step": "call_normal",
                        "type": "call",
                        "step-inputs": {
                            "indexedbam": {
                                "type": "bam",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "alignedbamfile",
                                            "type": "bam",
                                            "step": "align_normal",
                                            "tag": "normal"
                                        }
                                    }
                                }
                            },
                            "reference": {
                                "type": "reference",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "ref",
                                            "type": "reference",
                                            "step": "input-step",
                                            "tag": "input"
                                        }
                                    }
                                }
                            }
                        },
                        "step-outputs": {
                            "calls": {
                                "type": "vcf"
                            }
                        }
                    }
                }
            },
            "untagged": {
                "steps": {
                    "1": {
                        "step": "detect_mutation",
                        "type": "joint_call",
                        "step-inputs": {
                            "normal_tag": {
                                "type": "vcf",
                                "mapping": {
                                    "provided": "#normal",
                                    "candidates": {
                                        "1": {
                                            "id": "calls",
                                            "type": "vcf",
                                            "step": "call_normal",
                                            "tag": "normal"
                                        },
                                        "2": {
                                            "id": "calls",
                                            "type": "vcf",
                                            "step": "call_normal",
                                            "tag": "normal"
                                        }
                                    }
                                }
                            },
                            "tumour_tag": {
                                "type": "vcf",
                                "mapping": {
                                    "provided": "#tumour",
                                    "candidates": {
                                        "1": {
                                            "id": "calls",
                                            "type": "vcf",
                                            "step": "call_tumour",
                                            "tag": "tumour"
                                        },
                                        "2": {
                                            "id": "calls",
                                            "type": "vcf",
                                            "step": "call_tumour",
                                            "tag": "tumour"
                                        }
                                    }
                                }
                            },
                            "references": {
                                "type": "reference",
                                "mapping": {
                                    "provided": "",
                                    "candidates": {
                                        "1": {
                                            "id": "ref",
                                            "type": "reference",
                                            "step": "input-step",
                                            "tag": "input"
                                        }
                                    }
                                }
                            }
                        },
                        "step-outputs": {
                            "joint calls": {
                                "type": "vcf"
                            }
                        }
                    }
                }
            }
        }
    }
}""")


class TumourNormalPipeline(unittest.TestCase):

  def test_graph(self):
    translator = PipelineTranslator(debug=False)
    translator.translate_string(_yml)
    translation = translator.pipeline()
    # print('/\\'*40)
    # print(translation)
    # print('/\\'*40)
    # self.assertTrue(True)
    tr_json = json.loads(translation)
    self.assertTrue(tr_json == _expected)


if __name__ == '__main__':
    unittest.main()


