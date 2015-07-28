from interface import execute, Result
from adaptive_executor import adaptive_execute
from element_types import displacement, traction, free_slip, crack,\
    slip, mixed
from meshing import line, circle, sphere

import logging
logging.getLogger('elastic').addHandler(logging.NullHandler())
