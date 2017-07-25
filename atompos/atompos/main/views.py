import json as simplejson
from django.http import HttpResponse, Http404
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from atompos.main import settings
from util import get_atom_data, get_atom_data_atb

def index(request):
  return render(request, 'index.html')

def _get_positions(request, position_function):
  if request.method != 'POST':
    raise Http404

  params = request.POST.dict()
  if 'csrfmiddlewaretoken' in params:
    params.pop('csrfmiddlewaretoken')

  pos = position_function(params)
  if "error" in pos:
    res = pos
    res.update({"version": settings.VERSION})
  else:
    res = {
      "molecule": pos,
      "version": settings.VERSION
    }
  return HttpResponse(
    simplejson.dumps(res, indent=2),
    content_type="application/json"
  )

@csrf_exempt
def generate(request):
  return _get_positions(request, get_atom_data)

@csrf_exempt
def load_atb(request):
  return _get_positions(request, get_atom_data_atb)
