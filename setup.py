import shutil, os, subprocess, sys, shutil

from pathlib import Path
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

EXTENSION_NAME = "biosouppy"

class CMakeExtension(Extension):
  def __init__(self, name, sourcedir=""):
    Extension.__init__(self, name, sources=[])
    self.sourcedir = Path(sourcedir).absolute()

class CMakeBuild(build_ext):
  def build_extension(self, ext):
    install_dir = Path(self.get_ext_fullpath(ext.name)).absolute().parent
    source_dir = Path(__file__).absolute().parent
    build_dir = Path(self.build_temp).absolute()

    if build_dir.exists():
      shutil.rmtree(build_dir)
      build_dir.mkdir()

    cfg_args = [
      f'-H {source_dir}',
      f'-B {build_dir}',

      f'-Dbiosouppy=1',
      f'-Dbiosouppy_ext_dir={install_dir}',

      f'-DCMAKE_BUILD_TYPE=Release',
    ]

    try:
      import ninja
      cfg_args += ["-GNinja"]
    except ImportError:
      pass 

    build_args = ['--target', 'install']
    if 'CMAKE_BUILD_PARALLEL_LEVEL' not in os.environ:
      if hasattr(self, "parallel") and self.parallel:
        build_args += [f'-j {self.parallel}']  

    if not build_dir.exists():
      os.makedirs(build_dir)

    subprocess.check_call(
      ['cmake', *cfg_args], cwd=build_dir, stderr=sys.stderr)

    subprocess.check_call(
      ['cmake', '--build', '.', *build_args], cwd=build_dir, stderr=sys.stderr) 

setup(
  name='biosouppy',
  version='0.10.0',
  description='fasta/fastq container type with compression',
  long_description='README.md',
  author='Robert Vaser',
  author_email='robert.vaser@gmail.com',
  maintainer='Tvrtko Brekalo',
  maintainer_email='brekalo.tvrtko@gmail.com',
  ext_modules=[CMakeExtension(EXTENSION_NAME)],
  cmdclass={'build_ext': CMakeBuild},
  zip_safe=False,
  python_requires=">=3.6",
  license="MIT"
)
