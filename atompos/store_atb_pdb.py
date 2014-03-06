#!/usr/bin/env python
import os
import sys

if __name__ == "__main__":
  os.environ.setdefault("DJANGO_SETTINGS_MODULE", "atompos.settings")

  from atompos.main.util import generate_atb_pdb, load_atb_pdb

  for atb_id in sys.argv[1:]:
    try:
      load_atb_pdb(atb_id, True)
    except:
      generate_atb_pdb(atb_id)
      load_atb_pdb(atb_id, True)
