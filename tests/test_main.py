""" Tests of wc_test command line interface (wc_test.__main__)

:Author: Name <email>
:Date: 2018-5-10
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_test import __main__
import capturer
import mock
import unittest


class TestCore(unittest.TestCase):

    def test_cli(self):
        with mock.patch('sys.argv', ['wc_test', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegexpMatches(context.Exception, 'usage: wc_test')

    def test_help(self):
        with __main__.App(argv=['--help']) as app:
            app.run()

    def test_command_1(self):
        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['command-1']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), '...')
                self.assertEqual(captured.stderr.get_text(), '...')

    def test_command_3(self):
        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['command-3',
                                    'arg-1 value',
                                    'arg-2 value',
                                    '--opt-arg-3', 'opt-arg-3 value',
                                    '--opt-arg-4', 'opt-arg-4 value']) as app:
                # run app
                app.run()

                # test that the arguments to the CLI were correctly parsed
                self.assertTrue(app.pargs.arg_1)
                self.assertTrue(app.pargs.arg_2)
                self.assertTrue(app.pargs.opt_arg_3)
                self.assertTrue(app.pargs.opt_arg_4)

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), '...')
                self.assertEqual(captured.stderr.get_text(), '...')