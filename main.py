#!/usr/bin/env python3
# Copyright 2026 Antoine Alexandre, Hospices Civils de Lyon (HCL)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys
import os

# Ensure the repo root is on sys.path so `lv_explorer` is importable
# whether the script is run directly or bundled by PyInstaller.
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from lv_explorer.run import main

if __name__ == "__main__":
    main()
