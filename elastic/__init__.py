from interface import execute, Result
from element_types import displacement, traction, free_slip, crack, slip
from meshing import line, circle, sphere

import logging
logging.getLogger('elastic').addHandler(logging.NullHandler())
