from pathlib import Path
import shutil
import pdb
import os

from setuptools import find_packages, setup
from setuptools.command.build_py import build_py as _build_py

ROOT = Path(__file__).resolve().parent
PACKAGE_NAME = "unduwave"
EXTERNAL_ROOT = ROOT / PACKAGE_NAME / "External-Software"

WHITELIST = [
	("WAVE", "bin"),
	("WAVE", "stage"),
	("UNDUMAG", "bin"),
	("UNDUMAG", "stage"),
]

EXTRA_ROOT = ROOT / PACKAGE_NAME 
EXTRA_FOLDERS = [
	"MATERIAL_FILES",
	"UNDWAVE_IN_FILES",
]

class build_py(_build_py):
	def run(self):
		super().run()
		self.copy_external_files()
		self.copy_extra_folders()

	def copy_extra_folders(self):
		for rel_dir in EXTRA_FOLDERS:
			src_dir = EXTRA_ROOT / rel_dir
			dst_dir = Path(self.build_lib) / "unduwave" / rel_dir
			shutil.copytree(src_dir, dst_dir, dirs_exist_ok=True)

	def copy_external_files(self):
		for software, subdir in WHITELIST:
			src_dir = EXTERNAL_ROOT / software / subdir

			if not src_dir.exists():
				continue

			if subdir == 'stage' : 
				subdir='stage_tmp'

			dst_dir = (
				Path(self.build_lib)
				/  PACKAGE_NAME
				/ "External-Software"
				/ software
				/ subdir
			)

			for src_path in src_dir.rglob("*"):
				if src_path.is_file():
					rel = src_path.relative_to(src_dir)
					dst_path = dst_dir / rel
					dst_path.parent.mkdir(parents=True, exist_ok=True)
					shutil.copy2(src_path, dst_path)


setup(
	name="unduwave",
	version="0.8.0",
	package_dir={"": "."},
	packages=find_packages(where=".", include=["unduwave*"]),
	include_package_data=False,
	cmdclass={"build_py": build_py},
)