# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 27/09/2021
#
# Copyright : see License.txt
#
# Non-regression tests

import os
import subprocess

subprocess.run(["./../../../INSTALL-DIR/bin/Z_texture"], check=True)
subprocess.run(["./../../../INSTALL-DIR/bin/Main_pour_test"], check=True)
subprocess.run(["./../../../INSTALL-DIR/bin/Performance"], check=True)
subprocess.run(["./../../../INSTALL-DIR/bin/Performance_1"], check=True)
