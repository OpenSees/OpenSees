#!/usr/bin/env python
from __future__ import print_function
import os, sys
import json
try:
  import requests
except:
  raise RuntimeError('Run "sudo easy_install requests" before running this script')


#
#  If requests does not exist use sudo easy_install requests to install it.
#
host = os.getenv('SAWS_HOST')
if not host:
  host = 'localhost'
port = os.getenv('SAWS_PORT')
if not port:
  port = '8080'
url = 'http://'+host+':'+port+'/SAWs'


r = requests.get(url)
j = json.loads(r.content)

#  Example that access the functions in the stack
j = j['directories']['SAWs_ROOT_DIRECTORY']['directories']['PETSc']['directories']['Stack']['variables']['functions']['data']

print(j)


