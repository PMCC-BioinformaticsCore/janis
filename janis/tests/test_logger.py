from unittest import TestCase
from janis.utils.logger import Logger, LogLevel, _bcolors


class TestLogLevel(TestCase):
    def test_color_get_critical(self):
        self.assertEqual(_bcolors.FAIL, LogLevel.get_color(LogLevel.CRITICAL))

    def test_color_get_warning(self):
        self.assertEqual(_bcolors.WARNING, LogLevel.get_color(LogLevel.WARNING))

    def test_color_get_info(self):
        self.assertEqual(_bcolors.OKBLUE, LogLevel.get_color(LogLevel.INFO))

    def test_color_get_debug(self):
        self.assertEqual(_bcolors.ENDC, LogLevel.get_color(LogLevel.DEBUG))

    def test_color_get_other(self):
        self.assertEqual(_bcolors.ENDC, LogLevel.get_color(0))

    def test_str_get_critical(self):
        self.assertEqual("CRITICAL", LogLevel.get_str(LogLevel.CRITICAL))

    def test_str_get_warning(self):
        self.assertEqual("WARN", LogLevel.get_str(LogLevel.WARNING))

    def test_str_get_info(self):
        self.assertEqual("INFO", LogLevel.get_str(LogLevel.INFO))

    def test_str_get_debug(self):
        self.assertEqual("DEBUG", LogLevel.get_str(LogLevel.DEBUG))

    def test_str_get_other(self):
        self.assertEqual("", LogLevel.get_str(0))


class TestLogger(TestCase):
    def test_mute_unmute(self):
        Logger.set_console_level(LogLevel.DEBUG)
        self.assertIsNotNone(Logger.CONSOLE_LEVEL)
        Logger.mute()
        self.assertEqual(Logger.CONSOLE_LEVEL, LogLevel.CRITICAL)
        Logger.unmute()
        self.assertEqual(Logger.CONSOLE_LEVEL, LogLevel.DEBUG)

    def test_set_console_level(self):
        Logger.set_console_level(LogLevel.DEBUG)
