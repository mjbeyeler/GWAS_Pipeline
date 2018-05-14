import logging
import pysnptools.test
import unittest


if __name__ == '__main__':

    logging.basicConfig(level=logging.WARN)

    suites = unittest.TestSuite([pysnptools.test.getTestSuite()])
    suites.debug

    r = unittest.TextTestRunner(failfast=False)
    r.run(suites)

    logging.info("done with testing")
