from os import environ
from os.path import dirname
from sys import path

path = dirname(__file__)
if path not in path:
    path.append(path)
print path

environ['DJANGO_SETTINGS_MODULE'] = 'atompos.settings'

from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
