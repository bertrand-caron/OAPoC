#!/usr/bin/env python
import os
import sys

if __name__ == "__main__":
  os.environ.setdefault("DJANGO_SETTINGS_MODULE", "atompos.settings")

  from atompos.main.util import load_atb_pdb

  for atb_id in sys.argv[1:]:
    load_atb_pdb(atb_id, True)
