import json as simplejson
from django.http import HttpResponse, Http404
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from atompos.main import settings
from util import get_atom_pos, get_positions_atb

def index(request):
  return render(request, 'index.html')

@csrf_exempt
def generate(request):
  if request.method != 'POST':
    raise Http404

  params = request.POST.dict()
  if 'csrfmiddlewaretoken' in params:
    params.pop('csrfmiddlewaretoken')

  pos = get_atom_pos(params)
  pos.update({'version': settings.VERSION})
  return HttpResponse(
    simplejson.dumps(pos, indent=2, default=(lambda o: o.__dict__)),
    mimetype="application/json"
  )

@csrf_exempt
def load_atb(request):
  if request.method != 'POST':
    raise Http404

  params = request.POST.dict()
  if 'csrfmiddlewaretoken' in params:
    params.pop('csrfmiddlewaretoken')

  pos = get_positions_atb(params)
  pos.update({'version': settings.VERSION})
  return HttpResponse(
    simplejson.dumps(pos, indent=2, default=(lambda o: o.__dict__)),
    mimetype="application/json"
  )

