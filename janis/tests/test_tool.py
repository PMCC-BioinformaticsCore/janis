import unittest

from janis import Logger


class TestCommandTool(unittest.TestCase):
    def setUp(self):
        Logger.mute()

        from janis.unix.tools.tar import Tar

        self.tarTool = Tar()

        Logger.unmute()
