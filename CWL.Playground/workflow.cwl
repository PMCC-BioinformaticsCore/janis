#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0

requirements:
 - class: InlineJavascriptRequirement
 - class: ScatterFeatureRequirement
 - class: StepInputExpressionRequirement
 - class: SubworkflowFeatureRequirement

inputs:
  message:
    type: string

outputs:

steps:
  echo:
    run: echotool.cwl
    in: message

