from time import time
from logging import getLogger

# From: https://djangosnippets.org/snippets/1866/
def sizify(value):
    """
    Simple kb/mb/gb size snippet
    """
    #value = ing(value)
    if value < 512:
        ext = 'B'
    elif value < 512000:
        value = value / 1024.0
        ext = 'kB'
    elif value < 4194304000:
        value = value / 1048576.0
        ext = 'MB'
    else:
        value = value / 1073741824.0
        ext = 'GB'
    return '%s %s' % (str(round(value, 2)), ext)

class LoggingMiddleware(object):
    def __init__(self):
        # arguably poor taste to use django's logger
        self.logger = getLogger('django')

    def process_request(self, request):
        request.timer = time()
        return None

    def process_response(self, request, response):
        if not hasattr(request, 'timer'):
            request.timer = time()

        self.logger.info(
            '%s %s %s %s [%s] (%.0f ms)',
            request.META["SERVER_PROTOCOL"],
            request.method,
            request.get_full_path(),
            response.status_code,
            sizify(len(response.content)),
            (time() - request.timer) * 1000.
        )
        return response