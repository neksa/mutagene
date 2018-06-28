
from tests import *


class TestExample(unittest.TestCase):

    def test_pass(self):
        self.assertEqual(True, True)

    @pytest.mark.xfail
    def test_should_fail(self):
        self.assertEqual(False, True)
