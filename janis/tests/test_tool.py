import unittest

from janis import Logger


class TestCommandTool(unittest.TestCase):

    def setUp(self):
        Logger.mute()

        from janis.unix.tools.tar import Tar
        self.tarTool = Tar()

        Logger.unmute()

    def test_with_docker(self):
        tool_output = self.tarTool.cwl(with_docker=True).get_dict()

        self.assertIn("requirements", tool_output)
        reqs = tool_output["requirements"]
        self.assertIn("DockerRequirement", reqs)
        docker = reqs["DockerRequirement"]
        self.assertIn("dockerPull", docker)
        self.assertEqual(docker["dockerPull"], "ubuntu:latest")

    def test_without_docker(self):
        tool_output = self.tarTool.cwl(with_docker=False).get_dict()
        self.assertIn("requirements", tool_output)
        reqs = tool_output["requirements"]
        self.assertNotIn("DockerRequirement", reqs)
