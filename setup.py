from setuptools import setup
import setuptools.command.test

class TestCommand(setuptools.command.test.test):
  """ Setuptools test command explicitly using test discovery. """

  def _test_args(self):
    yield 'discover'
    for arg in super(TestCommand, self)._test_args():
      yield arg


setup(name='mouseAnno',
      version='0.1',
      description='Tools for parsing and testing gff annotation files',
      author='Gord Brown',
      author_email='gordon.brown@cruk.cam.ac.uk',
      license='MIT',
      packages=['mouse'],
      scripts=['bin/mouseAnno'],
      test_suite="tests",
      cmdclass={'test': TestCommand},
      zip_safe=True)
