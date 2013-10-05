import os
import sys

path = '/var/www/OAPoC'
if path not in sys.path:
    sys.path.append(path)

os.environ['DJANGO_SETTINGS_MODULE'] = 'atompos.settings'

import django.core.handlers.wsgi
application = django.core.handlers.wsgi.WSGIHandler()
