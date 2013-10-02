import simplejson
from django.http import HttpResponse, Http404
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from util import get_atom_pos

@csrf_exempt
def index(request):
  if request.method != 'POST':
    raise Http404

  params = request.POST.dict()
  try:
    params.pop('csrfmiddlewaretoken')
  except:
    pass

  pos = get_atom_pos(params)
  return HttpResponse(simplejson.dumps(pos, indent=2), mimetype="application/json")

def test(request):
  return render(request, 'test.html')

