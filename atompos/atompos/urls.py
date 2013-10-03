from django.conf.urls import patterns, include, url

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

urlpatterns = patterns('atompos.main.views',
    # Examples:
    # url(r'^$', 'atompos.views.home', name='home'),
    # url(r'^atompos/', include('atompos.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    # url(r'^admin/', include(admin.site.urls)),
    url(r'^OAPoC/$', 'index', name='index'),
    url(r'^OAPoC/test/', 'test', name='test'),
)
