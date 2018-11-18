"""
    Logger - Controls logging of the application
"""
from datetime import datetime

class _bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class LogLevel:
    NONE = 1
    CRITICAL = 2    # RED
    INFO = 3        # WHITE
    DEBUG = 4       # GREY

    @staticmethod
    def get_color(level: int):
        if level == LogLevel.CRITICAL:
            return _bcolors.FAIL

        if level == LogLevel.INFO:
            return _bcolors.OKBLUE

        if level == LogLevel.DEBUG:
            return _bcolors.ENDC

        return _bcolors.ENDC

    @staticmethod
    def get_str(level: int):
        if level == LogLevel.CRITICAL:
            return 'CRITICAL'

        if level == LogLevel.INFO:
            return 'INFO'

        if level == LogLevel.DEBUG:
            return 'DEBUG'

        return ''


class Logger:
    CONSOLE_LEVEL = LogLevel.CRITICAL
    DISK_LEVEL = LogLevel.DEBUG

    @staticmethod
    def log(message: str, level: int = LogLevel.DEBUG):
        if level == LogLevel.NONE:
            return

        if level <= Logger.CONSOLE_LEVEL:
            print(f"{LogLevel.get_color(level)} {Logger.get_prefix(level)}: {message}")


    @staticmethod
    def get_prefix(level: int):
        return f"{datetime.now().replace(microsecond=0).isoformat()} [{LogLevel.get_str(level)}]"
