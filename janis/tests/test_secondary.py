import unittest
from janis.translations.wdl import apply_secondary_file_format_to_filename


class TestApplySecondaryFileFormat(unittest.TestCase):

    def test_none(self):
        base = None
        sec = ".tsv"
        self.assertIsNone(apply_secondary_file_format_to_filename(base, sec))

    def test_append(self):
        base = "filename.ext"
        sec = ".tsv"
        self.assertEqual("filename.ext.tsv", apply_secondary_file_format_to_filename(base, sec))

    def test_one_extension(self):
        base = "filename.ext"
        sec = "^.tsv"
        self.assertEqual("filename.tsv", apply_secondary_file_format_to_filename(base, sec))

    def test_two_extensions(self):
        base = "filename.ext1.ext2"
        sec = "^^.tsv"
        self.assertEqual("filename.tsv", apply_secondary_file_format_to_filename(base, sec))

    def test_too_many_extensions(self):
        base = "filename.ext"
        sec = "^^.tsv"
        self.assertEqual("filename.tsv", apply_secondary_file_format_to_filename(base, sec))

    def test_way_too_many_extensions(self):
        base = "filename.ext"
        sec = "^^^^^^^^.tsv"
        self.assertEqual("filename.tsv", apply_secondary_file_format_to_filename(base, sec))

