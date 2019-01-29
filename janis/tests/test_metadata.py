import unittest
from datetime import date

from janis.utils.metadata import Metadata


class TestMetadata(unittest.TestCase):

    def test_update(self):

        m = Metadata()
        m2 = m.update(creator="c1")
        self.assertEqual(m.creator, "c1")
        self.assertEqual(m2.creator, "c1")

    def test_serialize(self):
        m = Metadata(creator="Michael Franklin")
        d = m.get_dict({"calculateChecksum": "ofThisDictionary"})
        self.assertIn("creator", d)
        self.assertEqual("e7443fb68e7ecf065100f014b2a2c16d9e357752", d["checksum"])

        # Check for all attributes that are none, and check they're not in the output
        for k, v in vars(m).items():
            if v is not None: continue
            self.assertNotIn(k, d)

    def test_date_generated(self):
        m = Metadata()
        d = m.get_dict({})
        self.assertIn("dateGenerated", d)
        self.assertEqual(date.today().strftime("%Y-%m-%d"), d["dateGenerated"])
