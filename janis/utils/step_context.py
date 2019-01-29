# from typing import Dict
#
# import json
#
# # from pipeline_definition.types.input_type import InputType
# from Workflow.step import Step
#
#
# class StepContext:
#   def __init__(self, step):
#     self.__step = step
#     self.__branch_outputs_stack = list()
#     self.__dependency_contexts = list()
#
#   def inherit_context_of_branch(self, prev_ctx):
#     outputs_stack_from_prev_context = prev_ctx.output_stack_of_branch()
#     if not outputs_stack_from_prev_context:
#       return
#     self.__branch_outputs_stack.extend(outputs_stack_from_prev_context)
#
#   def output_stack_of_branch(self):
#     output_stack = list()
#
#     provides_dict = {self.__step.id(): self.__step.provides()}
#
#     output_dict = {self.__step.tag(): provides_dict}
#     output_stack.append(output_dict)
#
#     output_stack.extend(self.__branch_outputs_stack)
#
#     return output_stack
#
#   def add_dependency_context_from(self, step_ctx):
#     dependency_context = step_ctx.provides()
#     self.__dependency_contexts.append(dependency_context)
#     # Dependency search is not only limited to last step in referred branch
#     self.__dependency_contexts.extend(step_ctx.output_stack_of_branch())
#
#   def provides(self):
#     provides_dict = {self.__step.id(): self.__step.provides()}
#     output_dict = {self.__step.tag(): provides_dict}
#     return output_dict
#
#   def provides_for(self, step):
#     provides = self.provides()
#
#     if not provides:
#       return []
#
#     tag_specific = provides.get(step.tag())
#
#     if not tag_specific:
#       return []
#
#     step_specific = tag_specific.get(step.id())
#     if not step_specific:
#       return []
#
#     return step_specific
#
#   def to_json(self):
#     return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)
#
#   def print(self):
#     print(self.to_json())
#
#   def map_input(self, step_input):
#
#     # mapping provided?
#     provided_mapping = self.__step.provided_value_for_requirement(step_input[Step.STR_ID])
#     if provided_mapping:
#       doc = {'provided': provided_mapping}
#     else:
#       doc = {'provided': ""}
#
#     if provided_mapping:
#       dependency_spec = Step.dependency_spec_from(provided_mapping)
#       candidates = self.find_match_for_input_dependency(step_input, dependency_spec)
#     else:
#       candidates = self.find_match_for_input(step_input)
#
#     if not candidates:
#       candidates['ERROR'] = "Failed to find an input candidate."
#
#     doc['candidates'] = candidates
#     return doc
#
#   def map_input_for_translation(self, step_input: InputType) -> Dict:
#
#     # mapping provided?
#     provided_mapping = self.__step.provided_value_for_requirement(step_input)
#
#     if provided_mapping:
#       dependency_spec = Step.dependency_spec_from(provided_mapping)
#       candidates = self.find_match_for_input_dependency(step_input, dependency_spec)
#     else:
#       candidates = self.find_match_for_input(step_input)
#
#     if not candidates:
#       raise Exception(f"Failed to find an input candidate for step: {self.__step.id()}, input requirement for: " +
#                       f"{step_input.type_name()}")
#
#     return candidates
#
#   def find_match_for_input_dependency(self, step_input, dependency_spec):
#
#     input_id = step_input[Step.STR_ID]
#     input_type = step_input[Step.STR_TYPE]
#     input_step_tag = self.__step.tag()
#
#     d_tag = dependency_spec['tag']
#
#     matches = {}
#     pref = 1
#
#     # Find the matching dependency in the stack
#     for dependency in self.__dependency_contexts:
#       matched = False
#       tag, step = next(iter(dependency.items()))
#       step_name, outputs = next(iter(step.items()))
#
#       if tag == d_tag:
#         # Found the matching dependency, now lets map the output
#
#         # Pass 1 to see if we have exact name / type match
#         for o in outputs:
#           o_id = o[Step.STR_ID]
#           o_type = o[Step.STR_TYPE]
#
#           # Name and Type match is heighest priority - conclusive
#           name = input_id
#           if name == o_id and input_type == o_type:
#             matches[pref] = self.__match_doc_for(o, step_name, input_step_tag)
#             pref = pref + 1
#             matched = True
#             break
#
#         if matched:
#           break
#
#         # Pass two is tag, tag_str type convention and type match
#         for o in outputs:
#           o_id = o[Step.STR_ID]
#           o_type = o[Step.STR_TYPE]
#
#           name = input_step_tag
#           if (o_id == name or o_id.startswith(name + "_")) and input_type == o_type:
#             matches[pref] = self.__match_doc_for(o, step_name, tag)
#             pref = pref + 1
#             matched = True
#             break
#
#         if matched:
#           break
#
#         # Pass three is type match
#         for o in outputs:
#           o_type = o[Step.STR_TYPE]
#
#           if input_type == o_type:
#             matches[pref] = self.__match_doc_for(o, step_name, tag)
#             pref = pref + 1
#
#         # Pass four is name match
#         for o in outputs:
#           o_id = o[Step.STR_ID]
#           if input_id == o_id:
#             matches[pref] = self.__match_doc_for(o, step_name, tag)
#             pref = pref + 1
#             continue
#
#     if matches:
#       return matches
#
#     return self.find_match_for_input(step_input)
#
#   def find_match_for_input(self, step_input: InputType):
#     # Finding the right input is based on the following criteria:
#     # 1. The output type of the nearest previous step matches
#     # 2. tags match                                        <<== This bit is not implemented yet
#     # 3. Search is limited if an input_scope is specified. <<== This bit is not implemented yet
#
#     step_tag = self.__step.tag()
#
#     matches = {}
#     pref = 1
#
#     # For each step in the stack, look at its provided outputs
#     for priority_entry in self.__branch_outputs_stack:
#       ostep_tag, ostep = next(iter(priority_entry.items()))
#       ostep_name, outputs = next(iter(ostep.items()))
#
#       for o_id, o_type in outputs.items():
#         if step_input == o_type: # and ostep_tag == step_tag:
#           matches[pref] = self.__match_doc_for(ostep_name, ostep_tag, o_id)
#           pref = pref + 1
#           break
#
#     return matches
#
#   @staticmethod
#   def __match_doc_for(name, tag, id):
#     doc = dict()
#     doc['step'] = name
#     doc['tag'] = tag
#     doc['id'] = id
#     return doc
