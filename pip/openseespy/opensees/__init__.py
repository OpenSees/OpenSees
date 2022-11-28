import sys

# only work for 64 bit system
if sys.maxsize < 2**31:
    raise RuntimeError('64 bit system is required')

# platform dependent
if sys.platform.startswith('linux'):

    try:
        from openseespylinux.opensees import *

    except:
        raise RuntimeError('Failed to import openseespy on Linux.')

elif sys.platform.startswith('win'):

    try:
        from openseespywin.opensees import *

    except:
        raise RuntimeError('Failed to import openseespy on Windows.')

elif sys.platform.startswith('darwin'):

    try:
        from openseespymac.opensees import *

    except:
        raise RuntimeError('Failed to import openseespy on Mac.')


else:

    raise RuntimeError(sys.platform+' is not supported yet')
