from os.path import join
from django.conf.urls import include, url
import atompos.main.views as views

RELATIVE_URL = r'^oapoc'

def relative_url(url):
    return join(RELATIVE_URL, url) if bool(url) else RELATIVE_URL

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

urlpatterns = [
    # Examples:
    # url(r'^$', 'atompos.views.home', name='home'),
    # url(r'^atompos/', include('atompos.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    # url(r'^admin/', include(admin.site.urls)),
    url(relative_url(r'$'), views.index, name='index'),
    url(relative_url(r'generate/'), views.generate, name='generate'),
    url(relative_url(r'loadATB/'), views.load_atb, name='load_atb'),
]
