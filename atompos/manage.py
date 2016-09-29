#!/usr/bin/env python
import os
import sys

SUPPORTED_DJANGO_VERSION = (1, 10)

if __name__ == "__main__":
    from django import VERSION
    if not all([(expected == current) for (expected, current) in zip(SUPPORTED_DJANGO_VERSION, VERSION)]):
        raise Exception(
            'ERROR: Wrong Django version. Use django-{0}'.format(
                '.'.join(map(str, SUPPORTED_DJANGO_VERSION)),
            ),
        )

    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "atompos.settings")

    from django.core.management import execute_from_command_line

    execute_from_command_line(sys.argv)
