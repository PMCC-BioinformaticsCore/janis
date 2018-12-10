import unittest

from Pipeline import String, Array, File, Int


class FileSubclass(File):
    def name(self):
        return "test_receivefrom_subclass"


class Test_ReceiveFrom(unittest.TestCase):

    def test_str_str(self):
        s1 = String()
        s2 = String()
        self.assertTrue(s2.can_receive_from(s1))

    def test_str_optstr(self):
        s1 = String(optional=False)
        s2 = String(optional=True)
        self.assertTrue(s2.can_receive_from(s1))

    def test_optstr_str(self):
        s1 = String(optional=True)
        s2 = String(optional=False)
        self.assertFalse(s2.can_receive_from(s1))

    def test_optstr_optstr(self):
        s1 = String(optional=True)
        s2 = String(optional=True)
        self.assertTrue(s2.can_receive_from(s1))

    def test_int_str(self):
        i = Int()
        s = String()
        self.assertFalse(s.can_receive_from(i))

    def test_arraystr_arraystr(self):
        ar1 = Array(String())
        ar2 = Array(String())
        self.assertTrue(ar2.can_receive_from(ar1))

    def test_arraystr_arrayoptstr(self):
        ar1 = Array(String())
        ar2 = Array(String(optional=True))
        self.assertTrue(ar2.can_receive_from(ar1))

    def test_arrayoptstr_arraystr(self):
        ar1 = Array(String(optional=True))
        ar2 = Array(String(optional=False))
        self.assertFalse(ar2.can_receive_from(ar1))

    def test_arrayoptstr_arrayoptstr(self):
        ar1 = Array(String(optional=True))
        ar2 = Array(String(optional=True))
        self.assertTrue(ar2.can_receive_from(ar1))

    def test_inheritance_forward(self):
        f1 = FileSubclass()
        f2 = File()
        self.assertTrue(f2.can_receive_from(f1))

    def test_inheritance_backward(self):
        f1 = File()
        f2 = FileSubclass()
        self.assertFalse(f2.can_receive_from(f1))