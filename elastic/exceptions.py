import traceback
import logging
logger = logging.getLogger(__name__)

def log_exceptions(f):
    def _wrapper(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            upper = ''.join(traceback.format_list(traceback.extract_stack()[:-1]))
            lower = '\n'.join(traceback.format_exc().split('\n')[1:-1])
            logger.error('Exception:\n' + upper + lower)
            raise
    return _wrapper

