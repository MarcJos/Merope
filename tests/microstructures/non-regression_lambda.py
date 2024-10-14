# -*- coding:utf8 -*-
#
# Copyright : see License.txt
#
# Non-regression tests

import os
import time
import sys
import importlib
import subprocess

def go_to_dir(name_dir):
    time.sleep(1)
    os.chdir(name_dir)
    print(name_dir)
    time.sleep(1)

def execute_tests(name_dir):
    go_to_dir(name_dir)
    subprocess.run(["python3", name_dir + ".py"])
    importlib.import_module(name_dir + ".test_vtk")
    go_to_dir("../")

list_of_test_names = [
    "parallel_gaussian",
    "prescribedField",
    "combineGeometryAndField",
    "texture"]

for test_name in list_of_test_names:
    execute_tests(test_name)



