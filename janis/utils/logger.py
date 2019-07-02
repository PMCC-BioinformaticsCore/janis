"""
    Logger - Controls logging of the application
"""
from datetime import datetime
from typing import Optional, TextIO


class _bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


class LogLevel:
    CRITICAL = 2  # RED
    WARNING = 3  # YELLOW
    INFO = 4  # WHITE
    DEBUG = 5  # GREY

    @staticmethod
    def get_color(level: int):
        if level == LogLevel.CRITICAL:
            return _bcolors.FAIL

        if level == LogLevel.WARNING:
            return _bcolors.WARNING

        if level == LogLevel.INFO:
            return _bcolors.OKBLUE

        if level == LogLevel.DEBUG:
            return _bcolors.ENDC

        return _bcolors.ENDC

    @staticmethod
    def get_str(level: int):
        if level == LogLevel.CRITICAL:
            return "CRITICAL"

        if level == LogLevel.WARNING:
            return "WARN"

        if level == LogLevel.INFO:
            return "INFO"

        if level == LogLevel.DEBUG:
            return "DEBUG"

        return ""


class Logger:
    CONSOLE_LEVEL: Optional[int] = LogLevel.INFO
    __TEMP_CONSOLE_LEVEL: Optional[int] = None
    WRITE_LEVEL: Optional[int] = LogLevel.DEBUG

    WRITE_LOCATION: Optional[str] = None
    __WRITE_POINTER: Optional[TextIO] = None

    @staticmethod
    def set_console_level(level: Optional[int]):
        Logger.CONSOLE_LEVEL = level

    @staticmethod
    def mute():
        """
        Decreases log level to critical until unmute is called
        """
        Logger.__TEMP_CONSOLE_LEVEL = Logger.CONSOLE_LEVEL
        Logger.set_console_level(LogLevel.CRITICAL)

    @staticmethod
    def unmute():
        """
        unmute the console, if __TEMP is None, should not do anything
        :return:
        """
        if Logger.__TEMP_CONSOLE_LEVEL is not None:
            Logger.set_console_level(Logger.__TEMP_CONSOLE_LEVEL)

    @staticmethod
    def set_write_level(level: Optional[int]):
        Logger.WRITE_LEVEL = level

    @staticmethod
    def set_write_location(location: str):
        if Logger.__WRITE_POINTER is not None:
            Logger.__WRITE_POINTER.close()

        if location is None:
            Logger.WRITE_LOCATION = None
            Logger.__WRITE_POINTER = None
            return

        Logger.WRITE_LOCATION = location
        Logger.__WRITE_POINTER = open(location, "w+")

    @staticmethod
    def close_file():
        if Logger.__WRITE_POINTER is not None:
            Logger.__WRITE_POINTER.close()

    @staticmethod
    def log(message: str, level: int = LogLevel.DEBUG):
        if level is None:
            # This is a developer error, we should never try to log with no level, it's purely for
            return

        m = f"{Logger.get_prefix(level)}: {message}"
        if Logger.CONSOLE_LEVEL is not None and level <= Logger.CONSOLE_LEVEL:
            print(LogLevel.get_color(level) + m + _bcolors.ENDC)

        if (
            Logger.WRITE_LEVEL is not None
            and level <= Logger.WRITE_LEVEL
            and Logger.__WRITE_POINTER is not None
        ):
            Logger.__WRITE_POINTER.write(m + "\n")

    @staticmethod
    def info(message: str):
        Logger.log(message, LogLevel.INFO)

    @staticmethod
    def warn(message: str):
        Logger.log(message, LogLevel.WARNING)

    @staticmethod
    def critical(message: str):
        Logger.log(message, LogLevel.CRITICAL)

    @staticmethod
    def log_ex(ex: Exception):
        Logger.critical(str(ex))

    @staticmethod
    def get_prefix(level: int):
        return f"{datetime.now().replace(microsecond=0).isoformat()} [{LogLevel.get_str(level)}]"
