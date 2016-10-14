test_wsgi:
	/usr/local/python27/bin/uwsgi --ini oapoc_uwsgi.ini
.PHONY: test_wsgi

first_setup: 
	# First, install django if you don't already hve it. Requires: Python 2.7+
	# Then, install open-babel (Macs: `brew install open-babel`)
	# Then ...
	python manage.py migrate
	python manage.py runserver
.PHONY: first_setup
