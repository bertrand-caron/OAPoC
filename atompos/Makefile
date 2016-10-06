test_wsgi:
	/usr/local/python27/bin/uwsgi --ini oapoc_uwsgi.ini
.PHONY: test_wsgi

first_setup: 
	# First, install django if you don't already hve it. Requires: Python 2.7+
	# Then, install open-babel (Macs: `brew install open-babel`)
	# Get your sjango Version for future reference
	python -c "import django; print(django.get_version())"
	# Then ...
	pip install sqlite3
	python manage.py migrate
	python manage.py runserver
