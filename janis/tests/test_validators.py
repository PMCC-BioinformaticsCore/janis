import unittest

from janis.utils.validators import Validators


class TestValidators(unittest.TestCase):

    def test_valid_identifiers(self):
        self.assertTrue(Validators.validate_identifier("test_workflow"))

    def test_invalid_identifiers(self):
        self.assertFalse(Validators.validate_identifier("test-workflow"))
