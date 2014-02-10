#!/usr/bin/env python
import os
import sys

if __name__ == "__main__":
  os.environ.setdefault("DJANGO_SETTINGS_MODULE", "atompos.settings")

  from atompos.main.util import get_positions_atb

  for atb_id in sys.argv[1:]:
    get_positions_atb({"molid": atb_id})
