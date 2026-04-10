import pytest

pytest.importorskip("flask", reason="webapp tests require [web] extras")
