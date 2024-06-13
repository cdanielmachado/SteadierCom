#!/usr/bin/env python

"""Tests for `steadiercom` package."""


import unittest

from steadiercom.cli import main_run

class TestSteadiercom(unittest.TestCase):
    """Tests for `steadiercom` package."""

    def setUp(self):
        """Set up test fixtures, if any."""
        pass

    def tearDown(self):
        """Tear down test fixtures, if any."""
        pass

    def test1(self):
        """Test something."""
        df = main_run(
            models=['tests/data/*.xml'],
            media='M9',
            mediadb='tests/data/media_db.tsv',
            w_e=0.002,
            w_r=0.2,
        )

        assert df is not None and len(df) > 0